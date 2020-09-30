import cell2cell as c2c
import numpy as np
import pandas as pd

import scipy
import json
import uuid
import os

# SETUP
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)  +  '/'
currentdir +=  '/'

data_folder = parentdir + '/Data/'

files = dict()
files['rnaseq'] = data_folder + '/RNA-Seq/Celegans_RNASeqData_Cell.xlsx'
files['metadata'] = data_folder + '/RNA-Seq/Celegans_cell_metadata.tsv'
files['ppi'] = data_folder + '/PPI-Networks/Celegans-Curated-LR-pairs.xlsx'
files['phenotype'] = data_folder + '/Digital-3D-Map/Celegans_Physical_Distances_Min.csv'

rnaseq_setup = dict()
rnaseq_setup['gene_col'] = 'gene_id'
rnaseq_setup['drop_nangenes'] = True
rnaseq_setup['log_transform'] = False

meta_setup = dict()
meta_setup['sample_col'] = '#SampleID'
meta_setup['group_col'] = 'Groups'

ppi_setup = dict()
ppi_setup['protein_cols'] = ['Ligand_WB', 'Receptor_WB']

cutoff_setup = dict()
cutoff_setup['type'] = 'constant_value'
cutoff_setup['parameter'] = 10 # TPM

analysis_setup = dict()
analysis_setup['communication_score'] = 'expression_thresholding'
analysis_setup['cci_score'] = 'bray_curtis'
analysis_setup['cci_type'] = 'undirected'

# Load Files
rnaseq_data = c2c.io.load_rnaseq(rnaseq_file=files['rnaseq'],
                                 gene_column=rnaseq_setup['gene_col'],
                                 drop_nangenes=rnaseq_setup['drop_nangenes'],
                                 log_transformation=rnaseq_setup['log_transform'],
                                 format='auto')

meta = c2c.io.load_metadata(metadata_file=files['metadata'],
                            rnaseq_data=rnaseq_data,
                            sample_col=meta_setup['sample_col'],
                            format='auto')

ppi_data = c2c.io.load_ppi(ppi_file=files['ppi'],
                           interaction_columns=ppi_setup['protein_cols'],
                           rnaseq_genes=list(rnaseq_data.index),
                           format='auto')

physical_distance = pd.read_csv(files['phenotype'], index_col=0)

# Analysis
included_cells = sorted(list(set(rnaseq_data.columns) & set(physical_distance.columns)))
excluded_cells = sorted(list(set(rnaseq_data.columns) - set(included_cells)))


def run_genetic_algorithm(ppi_data,
                          physical_distance,
                          included_cells,
                          population_size=200,
                          generations=200,
                          inc_percentage=0.05,
                          runs=None,
                          verbose=False
                          ):
    '''

    inc_percetange : float, 0.05 by default
        Min percentage of the objective function to improve in the last iteration with respect to the previous
        iteration. If this is not met, the runs are stopped.

    runs : int, None by default
        Number of runs of the GA. If runs is None, it iterates until the difference with the previous iteration
        is lower than 'imp_percentage' or 100 runs are performed.
    '''
    from pyevolve import G1DBinaryString
    from pyevolve import GSimpleGA
    from pyevolve import Selectors
    from pyevolve import Mutators

    print('Starting a new run of the genetic algorithm')

    results = dict()

    ref_distance_vector = scipy.spatial.distance.squareform(
        physical_distance.loc[included_cells, included_cells].values)

    theta_ppi_data = ppi_data.copy()

    runs_ = 1
    if runs is None:
        pass
    elif runs < 1:
        raise ValueError('Minimal number of runs is 1')

    old_obj = 0.1
    new_obj = old_obj * (1.1 + inc_percentage)
    while True:
        if runs is None:
            if (runs_ > 100) or ((new_obj - old_obj) / old_obj < inc_percentage):
                break
        else:
            if runs_ > runs:
                break
        # Sample only selected PPIs
        theta_ppi_data = theta_ppi_data.loc[theta_ppi_data.score == 1]

        bi_ppi_data = c2c.preprocessing.bidirectional_ppi_for_cci(ppi_data=theta_ppi_data, verbose=verbose)

        interaction_space = c2c.core.InteractionSpace(rnaseq_data=rnaseq_data[included_cells],
                                                      ppi_data=bi_ppi_data,
                                                      gene_cutoffs=cutoff_setup,
                                                      communication_score=analysis_setup['communication_score'],
                                                      cci_score=analysis_setup['cci_score'],
                                                      cci_type=analysis_setup['cci_type'],
                                                      verbose=verbose)

        def optimal_ppi_score(theta):
            # Inputs of this function
            theta_ppi_data['score'] = theta
            int_space = interaction_space  # InteractionSpace

            # Create bi_ppi_data to consider both directions (InteractionSpace previously created with bi_ppi_data)
            bi_ppi_data = c2c.preprocessing.bidirectional_ppi_for_cci(ppi_data=theta_ppi_data, verbose=verbose)
            int_space.ppi_data['score'] = bi_ppi_data['score']

            # Compute interactions
            int_space.compute_pairwise_cci_scores(use_ppi_score=True, verbose=verbose)

            # Distance matrix from analysis -> It considers only paracrine interaction (autocrine are excluded)
            result_matrix = interaction_space.distance_matrix.loc[included_cells, included_cells]

            # Compute correlations
            result_distance_vector = scipy.spatial.distance.squareform(result_matrix)
            corr = scipy.stats.spearmanr(result_distance_vector, ref_distance_vector)
            return abs(np.nan_to_num(corr[0]))

        # Genome instance
        genome = G1DBinaryString.G1DBinaryString(len(theta_ppi_data))

        # The evaluator function (objective function)
        genome.evaluator.set(optimal_ppi_score)
        genome.mutator.set(Mutators.G1DBinaryStringMutatorFlip)

        # Genetic Algorithm Instance
        ga = GSimpleGA.GSimpleGA(genome)
        ga.selector.set(Selectors.GTournamentSelector)
        ga.setPopulationSize(population_size)
        ga.setGenerations(generations)

        # Do the evolution, with stats dump
        # frequency of 10 generations
        ga.evolve(freq_stats=0)

        # Best result:
        best = ga.bestIndividual()
        theta_ppi_data['score'] = best

        # PRINT STATS:
        opt_score = optimal_ppi_score(best)
        drop_fraction = 1.0 - sum(best) / len(best)


        # Save PPI Networks
        tmp_ppi = ppi_data.assign(score=0)
        tmp_ppi.columns = ['A', 'B', 'score']
        tmp_ppi = tmp_ppi.merge(theta_ppi_data.loc[theta_ppi_data['score'] == 1],
                                how='outer',
                                on=['A', 'B'])

        tmp_ppi['score'] = tmp_ppi['score_y'].fillna(0)


        # Next iteration
        new_obj = opt_score
        if runs_ == 1:
            old_obj = new_obj / (1.1 + inc_percentage)
        else:
            old_obj = results['run{}'.format(runs_ - 1)]['obj_fn']

        # Save results
        results['run{}'.format(runs_)] = {'obj_fn': opt_score,
                                          'ppi_data': tmp_ppi['score'].values.tolist(),
                                          'drop_fraction': drop_fraction
                                          }

        runs_ += 1
    return results

if __name__ == "__main__":
    filename = data_folder + '/GA/' + str(uuid.uuid4())

    ga_results = run_genetic_algorithm(ppi_data,
                                       physical_distance,
                                       included_cells,
                                       population_size=200,
                                       generations=200,
                                       inc_percentage=0.025,
                                       runs=None,
                                       verbose=False
                                       )

    with open(filename + '.json', 'w') as fp:
        json.dump(ga_results, fp)

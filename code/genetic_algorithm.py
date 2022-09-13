import argparse

parser = argparse.ArgumentParser(description='Specify all inputs necessary for running the Genetic Algorithm Analysis.')
parser.add_argument("-s", "--score", dest="score", type=str, required=False, default="bray_curtis", choices=['bray_curtis', 'count', 'jaccard', 'icellnet'], help="CCI score to use")
parser.add_argument("-o", "--output_folder", dest="output_folder", type=str, required=False, default="GA", help="Full path to the folder where the results will be output (optional)")
parser.add_argument("-r", "--runs", dest="runs", type=int, required=False, default=100, help="Number of runs of the GA")
parser.add_argument("-c", "--cores", dest="cores", type=int, required=False, default=10, help="Number of cores to perform the all GA runs")

args = parser.parse_args()

import os

import cell2cell as c2c
import numpy as np
import pandas as pd

import scipy
import json
import uuid

# SETUP
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)  +  '/'
currentdir +=  '/'

data_folder = parentdir + '/data/'
output_folder = data_folder + '/' + args.output_folder + '/'
if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

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
analysis_setup['cci_score'] = args.score
analysis_setup['cci_type'] = 'undirected'

if args.score == 'icellnet':
    analysis_setup['communication_score'] = 'expression_product'

# Load Files
rnaseq_data = c2c.io.load_rnaseq(rnaseq_file=files['rnaseq'],
                                 gene_column=rnaseq_setup['gene_col'],
                                 drop_nangenes=rnaseq_setup['drop_nangenes'],
                                 log_transformation=rnaseq_setup['log_transform'],
                                 format='auto')
                                 
if args.score == 'icellnet':
    rnaseq_data = rnaseq_data.applymap(lambda x: np.log2(x+1))

meta = c2c.io.load_metadata(metadata_file=files['metadata'],
                            cell_labels=list(rnaseq_data.columns),
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
    
    
def run_ga(inputs):
    ga_results = run_genetic_algorithm(ppi_data=inputs['ppi_data'],
                                       physical_distance=inputs['physical_distance'],
                                       included_cells=inputs['included_cells'],
                                       population_size=inputs['population_size'],
                                       generations=inputs['generations'],
                                       inc_percentage=inputs['inc_percentage'],
                                       runs=inputs['runs'],
                                       verbose=inputs['verbose']
                                       )
    return ga_results

from multiprocessing import Pool

if __name__ == "__main__":
    
    inputs = {'ppi_data' : ppi_data, 'physical_distance' : physical_distance, 'included_cells' : included_cells,
              'population_size' : 200, 'generations' : 200, 'inc_percentage' : 0.025, 'runs' : None, 'verbose' : False}
    
    n_jobs = args.cores
    if n_jobs > args.runs:
        n_jobs = args.runs
    
    batches = int(np.ceil(args.runs/n_jobs))
    
    for i in range(batches):
        print(f'Running batch {i+1} of {batches}')
        with Pool(n_jobs) as p:
            input_ = [inputs]*n_jobs
            results = p.map(run_ga, input_)

        for r in results:
            filename = output_folder + str(uuid.uuid4())
            with open(filename + '.json', 'w') as fp:
                   json.dump(r, fp)
        print('')
        
            
    
    

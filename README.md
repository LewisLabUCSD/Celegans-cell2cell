# Evaluating cell-cell interactions and communication in C. elegans

## How to cite:

- Armingol E., Joshi C.J., Baghdassarian H., Shamie I., Ghaddar A., Chan J.,
 Her H.L., O’Rourke E.J., Lewis N.E. 
 [Inferring the spatial code of cell-cell interactions and communication across a whole animal body](https://doi.org/10.1101/2020.11.22.392217).
  *bioRxiv*, (2020). **DOI: 10.1101/2020.11.22.392217**

## Installation

All the analysis in this repository can be run on a **CodeOcean capsule** for reproducible results (estimated time: 2:30 hr): 

If you are interested in running everything locally, follow the instructions below.

### Installing Anaconda

Follow this tutorial to install [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation)

### Installing Github

Follow this tutorial to install [Github](https://gist.github.com/derhuerst/1b15ff4652a867391f03)

### Installing cell2cell
Create a new conda environment:
```
conda create -n cell2cell -y python=3.7 jupyter
```

Activate that environment:

```
conda activate cell2cell
```

Then, install all dependencies for reproducibility:
```
pip install 'numba==0.55.1'
pip install 'umap-learn==0.5.2'
pip install 'matplotlib==3.2.0'
pip install 'numpy==1.21.0'
pip install cell2cell
pip install git+https://github.com/BubaVV/Pyevolve
pip install tissue_enrichment_analysis
pip install matplotlib_venn
```

## Run/Explore the analyses

If the environment cell2cell is not active, activate it:

```
conda activate cell2cell
```

Then, for the respective analyses, a jupyter notebook is provided. Otherwise, instructions are detailed.

Jupyter notebooks can be open by executing the following command from the folder of this repository:

```
jupyter notebook
```

* Analyses have to be run in the order below (although we provided the results of each step, so analyses can be run skipping previous steps) :
1. [Generate list of ligands-receptors interactions from orthologs](./Notebooks/01.Generate-Celegans-LRs.ipynb).
Then, a manual-curation is needed. 
**This step can be skipped since we provided a [manual curated list](./Data/PPI-Networks/Celegans-Curated-LR-pairs.xlsx)**
2. [Compute intercellular distances and classify cell pairs by ranges of distances.](./Notebooks/02.Celegans-Cells-3D-Map.ipynb)
3. [Compute cell-cell interactions and communication for the curated ligand-receptor interactions.](./Notebooks/03.CCI-Curated-LRs.ipynb)
4. Run the genetic algorithm to select important ligand-receptor pairs for obtaining a better correlation
between CCI scores and intercellular distance:
    - From the main directory of this repository, run:
        ```
        python ./Notebooks/genetic_algorithm.py -s bray_curtis -o GA -r 100 -c 10
        ``` 
     *Note: This step can take between 1-2 days, depending on the number of iterations assigned
      in the nested for loops in the .py file. By default, it runs 100 times the GA, distributed in 10 cores (change -c 10 for another number of cores).
      
      **This step can be skipped since we provided [the results of 100 runs of the genetic algorithm using the Bray-Curtis-like score](./Data/GA/)**
5. [Examine results from the genetic algorithm, select ligand-receptor pairs and perform enrichment analysis on their functions.](./Notebooks/05.Examine-GA-Results.ipynb)
6. [Compute cell-cell interactions and communication for the GA-selected ligand-receptor interactions.](./Notebooks/06.CCI-Selected-LRs.ipynb)
7. [Perform permutation analyses on GA-selected ligand-receptor pairs.](./Notebooks/07.Permutation-Analysis-LRs.ipynb)
8. [Evaluate active ligand-receptor interactions along the body of C. elegans. Assess enrichment/depletion.](./Notebooks/08.Anteroposterior-Enrichment.ipynb)
9. [Evaluate enrichment/depletion of ligand-receptor pairs given their use in different ranges of distance.](./Notebooks/09.Distance-Ranges-Enrichment.ipynb)
10. [Evaluate enrichment of phenotypes on the genes in the GA-selected list of ligand-receptor pairs.](./Notebooks/10.Phenotype-Enrichment.ipynb)
11. [Generate UMAP plots based on Jaccard distance of pairs of cells given active LR pairs.](./Notebooks/11.CCC-UMAP-Visualization.ipynb)
12. Run a similar analysis to the one in step 4, but this time using the Count value as the CCI score: 
    - From the main directory of this repository, run:
        ```
        python ./Notebooks/genetic_algorithm.py -s count -o GA2 -r 100 -c 10
        ``` 
     *Note: This step can take between 1-2 days, depending on the number of iterations assigned
      in the nested for loops in the .py file. By default, it runs 100 times the GA, distributed in 10 cores (change -c 10 for another number of cores).
      
      **This step can be skipped since we provided [the results of 100 runs of the genetic algorithm using the LR count score](./Data/GA2/)**
13. [Examine results from the genetic algorithm (LR count as CCI score), select ligand-receptor pairs and perform enrichment analysis on their functions.](./Notebooks/12.Examine-GA2-Results.ipynb)
14. [Compare GA-based selection of LR pairs by using either Bray-Curtis-like score or LR count score](./Notebooks/14.CCI-Score-Comparison.ipynb)


**Disclaimer:** Results from the jupyter notebooks may change depending on the installed versions of the dependencies of the respective analyses.
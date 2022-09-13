# Inferring a spatial code of cell-cell interactions across a whole animal body

## How to cite:

- Armingol E., Ghaddar A., Joshi C.J., Baghdassarian H., Shamie I., Chan J.,
 Her H.L., O’Rourke E.J., Lewis N.E. 
 [Inferring a spatial code of cell-cell interactions across a whole animal bodys](https://doi.org/10.1101/2020.11.22.392217).
  *bioRxiv*, (2020). **DOI: https://doi.org/10.1101/2020.11.22.392217**

## Installation

All the analyses in this repository can be run on a **CodeOcean capsule** for reproducible results (estimated running time: 4:00 hr): https://doi.org/10.24433/CO.4688840.v2

If you are interested in running everything ***locally***, follow the ***instructions below***.

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

Then, install all dependencies:
```
pip install numba
pip install umap-learn
pip install 'matplotlib==3.2.0'
pip install 'cell2cell==0.6.0'
pip install git+https://github.com/BubaVV/Pyevolve
pip install tissue_enrichment_analysis
pip install matplotlib_venn
pip install 'xgboost==1.6.2'
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

### Main analyses
* Analyses have to be run in the order below (although we provided the results of each step, so analyses can be run skipping previous steps) :
1. [Generate list of ligands-receptors interactions from orthologs](./code/01.Generate-Celegans-LRs.ipynb).
Then, a manual-curation is needed. 
**This step can be skipped since we provided a [manual curated list](./data/PPI-Networks/Celegans-Curated-LR-pairs.xlsx)**
2. [Compute intercellular distances and classify cell pairs by ranges of distances.](./code/02.Celegans-Cells-3D-Map.ipynb)
3. [Compute cell-cell interactions and communication for the curated ligand-receptor interactions.](./code/03.CCI-Curated-LRs.ipynb)
4. Run the genetic algorithm to select important ligand-receptor pairs for obtaining a better correlation
between CCI scores and intercellular distance:
    - From the main directory of this repository, run:
        ```
        python ./code/genetic_algorithm.py -s bray_curtis -o GA-Bray-Curtis -r 100 -c 10
        ``` 
     *Note: This step can take between 1-2 days, depending on the number of iterations assigned
      in the nested for loops in the .py file. By default, it runs 100 times the GA, distributed in 10 cores (change -c 10 for another number of cores).
      
      **This step can be skipped since we provided [the results of 100 runs of the genetic algorithm using the Bray-Curtis score](./data/GA-Bray-Curtis/)**
5. [Examine results from the genetic algorithm, select ligand-receptor pairs and perform enrichment analysis on their functions.](./code/05.Examine-GA-Bray-Curtis.ipynb)
6. [Compute cell-cell interactions and communication for the GA-selected ligand-receptor interactions.](./code/06.CCI-Selected-LRs.ipynb)
7. [Perform permutation analyses on GA-selected ligand-receptor pairs.](./code/07.Permutation-Analysis-LRs.ipynb)
8. [Evaluate active ligand-receptor interactions along the body of C. elegans. Assess enrichment/depletion.](./code/08.Anteroposterior-Enrichment.ipynb)
9. [Evaluate enrichment/depletion of ligand-receptor pairs given their use in different ranges of distance.](./code/09.Distance-Ranges-Enrichment.ipynb)
10. [Evaluate enrichment of phenotypes on the genes in the GA-selected list of ligand-receptor pairs.](./code/10.Phenotype-Enrichment.ipynb)
11. [Generate UMAP plots based on Jaccard distance of pairs of cells given active LR pairs.](./code/11.CCC-UMAP-Visualization.ipynb)
12. Run a similar analysis to the one in step 4, but this time using the LR Count score as the CCI score: 
    - From the main directory of this repository, run:
        ```
        python ./Notebooks/genetic_algorithm.py -s count -o GA-LR-Count -r 100 -c 10
        ``` 
     *Note: This step can take between 1-2 days, depending on the number of iterations assigned
      in the nested for loops in the .py file. By default, it runs 100 times the GA, distributed in 10 cores (change -c 10 for another number of cores).
      
      **This step can be skipped since we provided [the results of 100 runs of the genetic algorithm using the LR count score](./data/GA-LR-Count/)**
13. [Examine results from the genetic algorithm (LR count as CCI score), select ligand-receptor pairs and perform enrichment analysis on their functions.](./code/13.Examine-GA-LR-Count.ipynb)
14. Run a similar analysis to the one in step 4, but this time using the ICELLNET score as the CCI score: 
    - From the main directory of this repository, run:
        ```
        python ./Notebooks/genetic_algorithm.py -s icellnet -o GA-ICELLNET -r 100 -c 10
        ``` 
     *Note: This step can take between 1-2 days, depending on the number of iterations assigned
      in the nested for loops in the .py file. By default, it runs 100 times the GA, distributed in 10 cores (change -c 10 for another number of cores).
      
      **This step can be skipped since we provided [the results of 100 runs of the genetic algorithm using the LR count score](./data/GA-ICELLNET/)**
15. [Examine results from the genetic algorithm (ICELLNET as CCI score), select ligand-receptor pairs and perform enrichment analysis on their functions.](./code/15.Examine-GA-ICELLNET.ipynb)
16. [Compare GA-based selection of LR pairs by using Bray-Curtis score, LR Count score, or ICELLNET score](./code/16.CCI-Score-Comparison.ipynb)
17. [Analyze spatial properties associated to the location type of each LR pair](./code/17.LR-Location-vs-CC-Distance-Heatmap.ipynb)

### Benchmarking analyses
1. [Generate CCI scores from Bray-Curtis, LR Count, and Smillie scoring functions](./code/benchmarking/01.Compute-Binary-based-scores.ipynb)
2. [Generate CCI scores from ICELLNET scoring function](./code/benchmarking/02.Compute-ICELLNET-scores.ipynb)
3. [Generate CCI scores from CellChat scoring function](./code/benchmarking/03.Compute-CellChat-scores.ipynb)
4. [Benchmarking of threshold values for binary-based methods](./code/benchmarking/04.Threshold-Benchmarking.ipynb)
5. [Benchmarking of all CCI-scores - Classifiers for distinguishing distance range between cells](./code/benchmarking/05.CCI-score-Benchmarking.ipynb)

**Disclaimer:** Figures from the jupyter notebooks may differ from those in the paper,
depending on the installed versions of the dependencies of the respective analyses.
The same might happen with certain results that depends on external tools.
**For ensuring the figures look the same, use the CodeOcean capsule instead.**
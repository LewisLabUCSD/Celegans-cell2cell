# Evaluating cell-cell interactions and communication in C. elegans

## How to cite:

- Armingol E., Joshi C.J., Baghdassarian H., Shamie I., Ghaddar A., Chan J.,
 Her H.L., O’Rourke E.J., Lewis N.E. 
 [Inferring the spatial code of cell-cell interactions and communication across a whole animal body](https://doi.org/10.1101/2020.11.22.392217).
  *bioRxiv*, (2020). **DOI: 10.1101/2020.11.22.392217**

## Installation
This tutorial works on Linux and macOS. All jupyter notebooks work on Windows, but the step 4 works only on Unix-based OS.

### Installing Anaconda

[Follow this tutorial to install Anaconda](https://docs.anaconda.com/anaconda/install/)

### Installing cell2cell
Create a new conda environment:
```
conda create -n cell2cell -y python=3.7.8 jupyter
```

Activate that environment:

```
conda activate cell2cell
```

Then, install cell2cell:
```
pip install 'cell2cell==0.3.1'
```

### Installing PyEvolve

*(This library is required in the script 4)*

[Install github](https://gist.github.com/derhuerst/1b15ff4652a867391f03)

Install a python 3 compatible version of PyEvolve:
```
pip install git+https://github.com/BubaVV/Pyevolve
```
### Installing Enrichment Analysis for C. elegans

*(This library is required in the notebook 10)*

Run:
```
pip install tissue_enrichment_analysis
```

### Installing umap

*(This library is required in the notebook 11)*

Run:
```
pip install umap-learn
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
        ./Notebooks/04.GA-LR-Selection.sh
        ``` 
     *Note: This step can take between 2-5 days, depending on the number of iterations assigned
      in the nested for loops in the .sh file. By default, it runs 100 times the GA, distributed in 10 cores.
      
      **This step can be skipped since we provided [the results of 100 runs of the genetic algorithm](./Data/GA/)**
5. [Examine results from the genetic algorithm, select ligand-receptor pairs and perform enrichment analysis on their functions.](./Notebooks/05.Examine-GA-Results.ipynb)
6. [Compute cell-cell interactions and communication for the GA-selected ligand-receptor interactions.](./Notebooks/06.CCI-Selected-LRs.ipynb)
7. [Perform permutation analyses on GA-selected ligand-receptor pairs.](./Notebooks/07.Permutation-Analysis-LRs.ipynb)
8. [Evaluate active ligand-receptor interactions along the body of C. elegans. Assess enrichment/depletion.](./Notebooks/08.Anteroposterior-Enrichment.ipynb)
9. [Evaluate enrichment/depletion of ligand-receptor pairs given their use in different ranges of distance.](./Notebooks/09.Distance-Ranges-Enrichment.ipynb)
10. [Evaluate enrichment of phenotypes on the genes in the GA-selected list of ligand-receptor pairs.](./Notebooks/10.Phenotype-Enrichment.ipynb)
11. [Generate UMAP plots based on Jaccard distance of pairs of cells given active LR pairs.](./Notebooks/11.CCC-UMAP-Visualization.ipynb)

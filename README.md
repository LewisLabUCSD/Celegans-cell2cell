# Evaluating cell-cell interactions and communication in C. elegans

## Installation
This tutorial works with Linux and macOS. For Windows, all jupyter notebooks work, the exception is the
step 4, with the genetic algorithm.

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
pip install 'cell2cell==0.2.2'
```

### Installing PyEvolve

[Install github](https://gist.github.com/derhuerst/1b15ff4652a867391f03)

Install a python 3 compatible version of PyEvolve:
```
pip install git+https://github.com/BubaVV/Pyevolve
```
### Installing obonet

Run:
```
pip install obonet
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

* Analyses have to be run in the order below:
1. [Generate list of ligands-receptors interactions from orthologs](./Notebooks/01.Generate-Celegans-LRs.ipynb).
Then, a manual-curation is needed. 
**This step can be omitted since we provided a [manual curated list](./Data/PPI-Networks/Celegans-Curated-LR-pairs.xlsx)**
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
5. [Examine results from the genetic algorithm, select ligand-receptor pairs and perform enrichment analysis on their functions.](./Notebooks/05.Examine-GA-Results.ipynb)
6. [Compute cell-cell interactions and communication for the GA-selected ligand-receptor interactions.](./Notebooks/06.CCI-Selected-LRs.ipynb)
7. [Perform permutation analyses on GA-selected ligand-receptor pairs.](./Notebooks/07.Permutation-Analysis-LRs.ipynb)
8. [Evaluate active ligand-receptor interactions along the body of C. elegans. Assess enrichment/depletion.](./Notebooks/08.Anteroposterior-Enrichment.ipynb)
9. [Evaluate enrichment/depletion of ligand-receptor pairs given their use in different ranges of distance.](./Notebooks/09.Distance-Ranges-Enrichment.ipynb)
10. [Evaluate enrichment of phenotypes on the genes in the GA-selected list of ligand-receptor pairs.](./Notebooks/10.Phenotype-Association.ipynb)

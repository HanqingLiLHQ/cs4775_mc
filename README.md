# Comparison of Computational Methods for Removing Redundancy in Motif Identification Process

This project provides Jupyter Notebooks and python files for performing hierarchical clustering and k-means clustering on various biological datasets, and runing (HMM) clustering, visualization of clusters, and distance calculations.

## Preparation

To run the modules and notebooks, ensure that you have Python(3.7 or later), Jupyter, biopython, matplotlib, numpy, scikit-learn (and maybe some other modules) installed. If not:

You can install python from their official website: 

`https://www.python.org/downloads/`

You can install Jupyter, biopython, matplotlib, scikit-learn (and some other modules) using pip:

```bash
pip install notebook
```
```bash
pip install biopython
```
```bash
pip install -U matplotlib
```
```bash
pip install numpy
```
```bash
pip install -U scikit-learn
```

Clone the repository to your local machine and get ready to run the project:

```bash
git clone https://github.com/HanqingLiLHQ/cs4775_mc
cd cs4775_mc
```

## Files on Main
### `FINAL_REPORT.pdf`
This is the final report for the project. It described the methods we used/developed for motif clustering, the results we got from these methods, and the discussion we had about each clustering method and its corresponding results. 

### `motif_distance.ipynb`

Run distance calculations using this notebook. It contains the necessary functions and code blocks to process the datasets and compute distances.

It was used to plot and compare different alignment methods and distance functions. Specific motifs were fetched from JASPAR, including similar ones, dissimilar ones, and gapped ones with different gap lengths. Among the functions implemented, "compare_distance" compares different position-wise distance calculation methods (used to plot figure 3 and 4 by comparing fetched motifs ppm1 and ppm2, ppm1 and ppm3, gapped_ppm1 and gapped_ppm2). "compare_align" compares the overlap and expand strategy for motif alignment, giving rise to figure 2 by comparing ppm1 and ppm2, ppm1 and ppm3.

### `distance_based_clustering.ipynb`

Use this notebook to run hierarchical and k-means clustering on the following datasets:
- Gapped vertebrates
- fungi_large
- fungi
- vertebrates_large

Each dataset has a dedicated code block. To switch between hierarchical and k-means clustering, comment out the code for the clustering algorithm you do not want to use, and uncomment the one you do.

### `hmm_clustering.py`

Run this python module to perform HMM clustering specifically on the `fungi` dataset.

```bash
python3 hmm_clustering.py
```

### `cluster_graph.py`

Run this python module for visualizing the clusters formed by the clustering algorithms. 
```bash
python3 cluster_graph.py
```
It will output the following folders:
- Motifs Graph: This folders contains motifs weblogos for JASPAR dataset fungi.

- Nature Clusters: This folder contains JASPAR's clusters on dataset fungi.

- Hierarchical Clusters: This folder contains graphing results for hierarchical clustering on dataset fungi.

- Kmeans Clusters: This folder contains graphing results for kmeans clustering on dataset fungi.

- HMM Clusters: This folder contains graphing results for hmm clustering on dataset fungi.

## Folders on Main

### `Dataset`

This folder contains all the datasets we used for this project.

### `Distance_Based_Clustering`

This folder contains all the methods we used for hierarchical and k-means clustering.

### `HMM_Clustering`

This folder contains all the methods we used for HMM Clustering

### `Evaluation`

This folder contains all the scoring methods and visualization tools we uesd for clustering results.

## Acknowledgements

Thanks to Hanqing Li, Esther Yu, Rainney Wan, Sicheng Ma, and Jingyu Xu for investing time in developing and analyzing these amazing clustering methods🥳🥳🥳🥳🥳

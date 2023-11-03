This folder contains the implementation of hierarchical agglomerative clustering methods. 

hierarchical_clustering.py
  Several linkage methods are build into biopython:
  1, Single-linkage clustering:
    The two branches expanding from the same node should have the same length while building the tree. 
    Thus all the initial elements are on the same plane (all rank 0). All the distances can be reconstructed based on rank.
    Step 1: from the current distance matrix, cluster the two elements with the shortest distance
    Step 2: update an internal node with rank = shortest_distance / 2, calculate the other distances based on rank accordingly. 
    Step 3: update the distance matrix by changing the distance to the new cluster as the minimum distance to the two elements. 

    it might be good to distinguish between closely related elements. 

    https://en.wikipedia.org/wiki/Single-linkage_clustering


  2, Complete-linkage clustering:
    The two branches expanding from the same node should have the same length while building the tree. 
    Thus all the initial elements are on the same plane (all rank 0). All the distances can be reconstructed based on rank.
    The only difference to single-linkage clustering is that it uses the maximum distance to update the distance matrix. 
    This would be good for avoiding the chaining of clusters: clusters formed from single-linkage clustering might be too close to each other. 
    Step 1: from the current distance matrix, cluster the two elements with the shortest distance
    Step 2: update an internal node with rank = shortest_distance / 2, calculate the other distances based on rank accordingly. 
    Step 3: update the distance matrix by changing the distance to the new cluster as the maximum distance to the two elements.  

    It might be good to cluster highly different clusters. 

    https://en.wikipedia.org/wiki/Complete-linkage_clustering

  3，Average-linkage clustering (UPGMA):
    The two branches expanding from the same node should have the same length while building the tree. 
    Thus all the initial elements are on the same plane (all rank 0). All the distances can be reconstructed based on rank.
    Step 1: from the current distance matrix, cluster the two elements with the shortest distance
    Step 2: update an internal node with rank = shortest_distance / 2, calculate the other distances based on rank accordingly. 
    Step 3: update the distance matrix by changing the distance to the new cluster as the (weighted) average distance to the two elements. 
    The weight of the two elements based on how much motifs each cluster element contain.     

    It it just changing the step to update the distance matrix. 
    It would be good for molecular clock assumption, giving the same weight to actually all the initial elements. 
    UPGMA would take into evolutionary aspect while complete and single does not. 

  4, The WPGMA method is simply changing the distance update to be the sum of the two distance divided by 2.
    It would not consider the size of the clusters but just consider them the same. 
  https://en.wikipedia.org/wiki/WPGMA



Above four tree-building methods are building trees with the same branch-length, and 
all the internal nodes have fixed ranks. Therefore, we just need to cut the tree at
a given rank, and clusters will be formed.
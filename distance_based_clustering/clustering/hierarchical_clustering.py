from Bio.Cluster import treecluster
import numpy as np
import os


def hierarchical_clustering(nclusters, dm, id_list):
    # Using Bio.Cluster.treecluster to do hierarchical clustering based on the distance matrix
    # possible parameters for method:
    # method = "s" : pairwise single-linkage clustering
    # method = "m" (default): pairwise complete-linkage clustering
    # method = "a": pairwise average-linkage clustering
    # method = "c": pairwise centroid-linkage clustering (cannot be calculated from the distance matrix)
    tree = treecluster(data=None, distancematrix=np.array(dm))

    # return the number of cluster each id belongs to
    cluster_id_list = tree.cut(nclusters)
    # assemble the clusters with motif_ids. Build up a dictionary where cluster ids are the keys,
    # and the lists containing motif_ids are values

    # initialize the clusters map
    clusters_map = {}
    for cluster_id in range(nclusters):
        clusters_map[cluster_id] = []

    # iterate through the cluster_id list to classify the motifs
    for i in range(len(cluster_id_list)):
        clusters_map[cluster_id_list[i]].append(id_list[i])

    clusters = list(clusters_map.values())

    return clusters

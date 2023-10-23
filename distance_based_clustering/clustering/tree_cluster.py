from Bio.Cluster import treecluster
import sys
sys.path.insert(1,'../distance_matrix_calculation')
import pearson_distance
import numpy as np
import os

# calculate the distance matrix and the ids
dm,id_list = pearson_distance.calculate_distance_matrix()

# build up the tree
tree = treecluster(data= None, distancematrix = np.array(dm))

# return the number of cluster each id belongs to
cluster_id_list = tree.cut(nclusters = 5)
# assemble the clusters with motif_ids. Build up a dictionary where [0-4] are the keys,
# and the lists containing motif_ids are values
clusters_map = {}
clusters_map[0] = []
clusters_map[1] = []
clusters_map[2] = []
clusters_map[3] = []
clusters_map[4] = []
for i in range(len(cluster_id_list)):
   clusters_map[cluster_id_list[i]].append(id_list[i])

clusters = list(clusters_map.values())

print(clusters)
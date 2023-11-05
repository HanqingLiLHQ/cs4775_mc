import sys
sys.path.insert(1,'../dataset')
import dbm
import dataset_yeast_core
from  clustering import hierarchical_clustering
from  clustering import kmeans, kmeans_self_defined_dist
from  distance_matrix_calculation import pearson_distance

dataset = dataset_yeast_core.DNABindingMotifs()
# dataset.mmm is a dictionary with format {"id": Bio.motifs.Motif object}
# where each motif matrix can be retrieved by motif.pssm
print(dataset.mmm)
print("Finish printing dataset \n")

# Hierarchical clustering
dm, idlist = pearson_distance.calculate_distance_matrix(dataset)
# # distance matrix
# print(dm)
# print("Finish printing distance matrix \n")
# # motif's order in the distance matrix
# print(idlist)
# print("Finish printing idlist \n")

# clusters = hierarchical_clustering.hierarchical_clustering(5, dm, idlist)
# # resulting cluters after using clustering method
# print(clusters)



# K-means clustering
num_clusters = 3  # Define the number of clusters
clusters = kmeans_self_defined_dist.kmeans_motifs(dm, idlist, num_clusters)
print(clusters)
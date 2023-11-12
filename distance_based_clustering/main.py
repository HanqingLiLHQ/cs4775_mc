import sys
sys.path.insert(1,'../dataset')
import dbm
import dataset_yeast_mini
from  clustering import hierarchical_clustering
from  clustering import kmeans, kmeans_self_defined_dist
from  distance_matrix_calculation import pearson_distance
from  distance_matrix_calculation import jensen_distance

dataset = dataset_yeast_mini.DNABindingMotifs()
# dataset.mmm is a dictionary with format {"id": Bio.motifs.Motif object}
# where each motif matrix can be retrieved by motif.pssm
print(dataset.mmm)
print("Finish printing dataset \n")
dm, idlist = pearson_distance.calculate_distance_matrix(dataset)
# dm, idlist = jensen_distance.calculate_distance_matrix(dataset)


# Hierarchical clustering
# clusters = hierarchical_clustering.hierarchical_clustering(5, dm, idlist)
# # resulting cluters after using clustering method
# print(clusters)


# K-means clustering
# One advantage of K-means is it can customize number of clusters needed
num_clusters = 3  # Define the number of clusters
clusters = kmeans_self_defined_dist.kmeans_motifs(dm, idlist, num_clusters)
print(clusters)
print("Finish printing result by K-means \n")
import sys
sys.path.insert(1,'../dataset')
import dbm
from  clustering import hierarchical_clustering
from  distance_matrix_calculation import pearson_distance

dataset =   dbm.DNABindingMotifs()
dm, idlist = pearson_distance.calculate_distance_matrix(dataset)
clusters = hierarchical_clustering.hierarchical_clustering(5, dm, idlist)

print(clusters)
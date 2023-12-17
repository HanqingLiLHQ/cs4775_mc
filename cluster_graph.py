import sys
import dbm
import os

dataset_path = os.path.join(os.getcwd(), "dataset")
sys.path.insert(0, dataset_path)
distance_based_clustering_path = os.path.join(os.getcwd(), "distance_based_clustering")
sys.path.insert(1, distance_based_clustering_path)
evaluate_path = os.path.join(os.getcwd(), "evaluate")
sys.path.insert(2, evaluate_path)
clustering_path = os.path.join(os.getcwd(), "distance_based_clustering", "clustering")
sys.path.insert(3, clustering_path)
import fungi
from column_distance import *
from distance_based_clustering.clustering import hierarchical_clustering
from distance_based_clustering.clustering import kmeans_self_defined_dist
from distance_based_clustering.distance_matrix_calculation import (
    pearson_distance,
)
from distance_based_clustering.distance_matrix_calculation import jensen_distance
from distance_based_clustering.distance_matrix_calculation import (
    dist_from_col,
)
from evaluate import ari_score
from evaluate import fmi_score
from evaluate import nmi_score
from evaluate import graph
from evaluate import motif_graph
from evaluate import cluster_graph
from evaluate import final_graph

################################################################################

# Helper Functions


def assign_unique_numbers_to_values(dataset):
    """
    Assign unique numbers to the values in the dataset dictionary.
    If two values are the same, they will be given the same number.

    :param dataset: dict, a dictionary with values that need unique numbers assigned.
    :return: dict, a dictionary with values replaced by unique numbers.
    """
    unique_values = set(dataset.values())  # Get the unique values
    value_to_number = {
        value: idx for idx, value in enumerate(unique_values)
    }  # Map each unique value to a number

    return {key: value_to_number[value] for key, value in dataset.items()}


def generate_labels_from_clusters(clusters):
    """
    Generate a dictionary of labels based on cluster membership.
    Items in the same list receive the same label, starting from 0.

    :param clusters: list of lists, where each sublist represents a cluster.
    :return: dict, a dictionary with items as keys and their respective labels as values.
    """
    label_dict = {}
    for label, cluster in enumerate(clusters):
        for item in cluster:
            label_dict[item] = label
    return label_dict


def create_cluster_dicts(clusters, dataset):
    """
    Create a list of dictionaries for each cluster.

    :param clusters: A list of clusters, each cluster being a list of item keys.
    :param dataset: The dataset containing the items and their corresponding data.
    :return: A list of dictionaries, each dictionary representing a cluster.
    """
    our_cluster = []
    for cluster in clusters:
        cluster_dict = {}
        for item in cluster:
            cluster_dict[item] = dataset.mmm[item]
        our_cluster.append(cluster_dict)
    return our_cluster


def most_similar_lists(A, B):
    def similarity(list1, list2):
        # Count common elements for similarity
        return len(set(list1) & set(list2))

    B_matched = []
    B_used = [False] * len(B)  # Keep track of which lists in B have been used

    for a in A:
        max_similarity = -1
        best_match_idx = -1

        for idx, b in enumerate(B):
            if not B_used[idx]:
                current_similarity = similarity(a, b)
                if current_similarity > max_similarity:
                    max_similarity = current_similarity
                    best_match_idx = idx

        B_matched.append(B[best_match_idx])
        B_used[best_match_idx] = True  # Mark this list as used

    return B_matched


def sort_sublists(A):
    for sublist in A:
        sublist.sort()
    return A


################################################################################

# Dataset: fungi

dataset = fungi.DNABindingMotifs()
# dataset.mmm is a dictionary with format {"id": Bio.motifs.Motif object}
print(dataset.mmm)
print("Finish printing fungi.mmm \n")

base_dir_motifs = "motifs_graphs"
motif_generator = motif_graph.MotifGraphGenerator(dataset.mmm, base_dir_motifs)
motif_generator.generate_motif_graphs()


dm, idlist = dist_from_col.calculate_distance_matrix(
    dataset, Euclidean_Distance, "expand", bg=[0.27, 0.23, 0.23, 0.27], average=np.mean
)

dataset_cc = dataset.cc
print(dataset_cc)
print("Finish printing dataset3 cc \n")

# base_dir3 = "nature_clusters3"
# generator3 = graph.MotifGraphGenerator(dataset3_cc, base_dir3)
# cluster_paths3 = generator3.generate_cluster_graphs()
# generator3.generate_final_composite_graph(cluster_paths3)

# Convert dataset3_cc into a list of lists, where each inner list contains the keys of the corresponding dictionary
dataset_cc_list = [[key for key in d] for d in dataset_cc]
dataset_cc_list.sort(key=len, reverse=True)
# dataset3_cc_list.reverse()
dataset_cc_list = sort_sublists(dataset_cc_list)
print(dataset_cc_list)
print("Finish printing dataset3 cc list \n")

# For clusters of motifs
base_dir_clusters = "nature_clusters3"
cluster_generator = cluster_graph.ClusterGraphGenerator(
    dataset_cc_list, base_dir_motifs, base_dir_clusters
)
cluster_paths = cluster_generator.generate_cluster_graphs()
final_graph.generate_final_composite_graph(base_dir_clusters, cluster_paths)


num_clusters = len(dataset_cc_list)  # Define the number of clusters
print(num_clusters)
print("Finish printing num_clusters3 \n")

# distance_update can be one of ["complete","single","UPGMA","WPGMA"]
cluster = hierarchical_clustering.hierarchical_clustering(
    num_clusters, dm, idlist, distance_update="complete"
)
cluster = sort_sublists(cluster)
print(cluster)
print("Finish printing result by hierarchical clustering \n")

cluster_match = most_similar_lists(dataset_cc_list, cluster)
print(cluster_match)
print("Finish printing result by matching \n")

base_dir_clusters_our = "predicted_clusters3"
cluster_generator_our = cluster_graph.ClusterGraphGenerator(
    cluster_match, base_dir_motifs, base_dir_clusters_our
)
cluster_paths_our = cluster_generator_our.generate_cluster_graphs()
final_graph.generate_final_composite_graph(base_dir_clusters_our, cluster_paths_our)

# cluster2 = kmeans_self_defined_dist.kmeans_motifs(dm2, idlist2, num_clusters2)
# print(cluster2)
# print("Finish printing result by K means \n")

# our_cluster3 = create_cluster_dicts(cluster3, dataset3)

# base_dir_our3 = "predicted_clusters3"
# generator_our3 = graph.MotifGraphGenerator(our_cluster3, base_dir_our3)
# cluster_paths_our3 = generator_our3.generate_cluster_graphs()
# generator_our3.generate_final_composite_graph(cluster_paths_our3)


dataset_dic = dict(sorted(generate_labels_from_clusters(dataset_cc_list).items()))
print(dataset_dic)
print("Finish printing true labels \n")

clusters_dic = dict(sorted(generate_labels_from_clusters(cluster).items()))
print(clusters_dic)
print("Finish printing predicted labels \n")

true_labels = list(dataset_dic.values())
predicted_labels = list(clusters_dic.values())

score_ari = ari_score.evaluate_clustering(true_labels, predicted_labels)
score_fmi = fmi_score.evaluate_clustering(true_labels, predicted_labels)
score_nmi = nmi_score.evaluate_clustering(true_labels, predicted_labels)

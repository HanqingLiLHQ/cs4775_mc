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
# base_dir_clusters = "nature_clusters3"
# cluster_generator = cluster_graph.ClusterGraphGenerator(
#     dataset_cc_list, base_dir_motifs, base_dir_clusters
# )
# cluster_paths = cluster_generator.generate_cluster_graphs()
# final_graph.generate_final_composite_graph(base_dir_clusters, cluster_paths)


num_clusters = len(dataset_cc_list)  # Define the number of clusters
print(num_clusters)
print("Finish printing num_clusters3 \n")

# distance_update can be one of ["complete","single","UPGMA","WPGMA"]
cluster = [
    ["MA0284.2"],
    ["MA0418.1"],
    ["MA0439.1"],
    ["MA0324.1"],
    ["MA0407.1"],
    ["MA0351.1"],
    ["MA0359.2"],
    ["MA0347.2"],
    ["MA0304.1"],
    ["MA0279.2"],
    ["MA0287.1"],
    ["MA0413.1"],
    ["MA0316.1"],
    ["MA0299.1"],
    ["MA0285.1"],
    ["MA0377.1"],
    ["MA0419.1"],
    ["MA0365.1"],
    ["MA0371.1"],
    ["MA0294.1"],
    ["MA0326.1"],
    ["MA0440.1"],
    ["MA0929.1"],
    ["MA0379.1", "MA0386.1"],
    ["MA0313.1", "MA0314.2"],
    ["MA0409.1", "MA0281.2"],
    ["MA0302.1", "MA0389.1", "MA0301.1"],
    ["MA1435.1", "MA0353.1"],
    ["MA0373.1", "MA0333.1", "MA0332.1"],
    ["MA0415.1", "MA0416.1"],
    ["MA0421.1", "MA0363.2"],
    ["MA0350.1", "MA0329.1", "MA0330.1"],
    ["MA0390.1", "MA0398.1"],
    ["MA0367.1", "MA0437.1"],
    ["MA0417.1", "MA0317.1", "MA0297.1", "MA0296.1"],
    ["MA0382.2", "MA0286.1"],
    ["MA0385.1", "MA0355.1"],
    ["MA0408.1", "MA0278.1"],
    ["MA0336.1", "MA0266.1"],
    ["MA0403.2", "MA0423.1", "MA0306.1", "MA0372.1"],
    ["MA0269.1", "MA0270.1"],
    ["MA0384.1", "MA0397.1", "MA0396.1"],
    ["MA0293.1", "MA1437.1"],
    ["MA0411.1", "MA0340.1", "MA0370.1"],
    [
        "MA0364.1",
        "MA1433.1",
        "MA1432.1",
        "MA0436.1",
        "MA0425.1",
        "MA0431.1",
        "MA0268.1",
        "MA0337.1",
    ],
    ["MA0414.1", "MA0400.1", "MA0282.1", "MA0391.1"],
    ["MA0368.1", "MA0381.1", "MA0361.1", "MA0429.1", "MA0320.1", "MA0412.2"],
    ["MA0288.1", "MA0303.2"],
    ["MA0420.1", "MA0352.2", "MA0360.1"],
    ["MA0434.1", "MA0435.1"],
    ["MA0395.1", "MA0291.1"],
    ["MA0323.1", "MA0362.1", "MA0399.1", "MA0349.1", "MA0404.1"],
    ["MA0426.1", "MA0433.1", "MA0356.1", "MA0369.1"],
    ["MA0276.1", "MA0331.1"],
    ["MA0402.2", "MA0394.1"],
    ["MA0310.1", "MA0357.1", "MA1436.1", "MA0376.1", "MA0305.1", "MA0406.1"],
    [
        "MA0388.1",
        "MA0277.1",
        "MA0327.1",
        "MA0334.1",
        "MA0343.1",
        "MA1434.1",
        "MA0295.1",
    ],
    ["MA1431.1", "MA0338.1", "MA0339.1", "MA0441.1", "MA0267.2"],
    ["MA0274.1", "MA0387.1"],
    [
        "MA0283.1",
        "MA0344.1",
        "MA0366.1",
        "MA0342.1",
        "MA0341.1",
        "MA0375.1",
        "MA0325.1",
        "MA0290.1",
        "MA0410.1",
    ],
    [
        "MA0428.1",
        "MA0311.1",
        "MA0392.1",
        "MA0430.1",
        "MA0422.1",
        "MA0358.1",
        "MA0374.1",
        "MA0401.1",
        "MA0424.1",
        "MA0280.1",
        "MA0275.1",
        "MA0292.1",
        "MA0308.1",
        "MA0380.1",
        "MA0405.1",
        "MA0312.2",
        "MA0354.1",
        "MA0438.1",
        "MA0432.1",
        "MA0348.1",
    ],
    ["MA0328.2", "MA0318.1", "MA0321.1", "MA0322.1"],
    ["MA0378.1", "MA0319.1", "MA0393.1"],
    ["MA0273.1", "MA0271.1", "MA0272.1"],
    ["MA0307.1", "MA0309.1", "MA0265.2", "MA0289.1", "MA0300.1"],
]
cluster = sort_sublists(cluster)
print(cluster)
print("Finish printing result by hierarchical clustering \n")

cluster_match = most_similar_lists(dataset_cc_list, cluster)
print(cluster_match)
print("Finish printing result by matching \n")

base_dir_clusters_our = "predicted_clusters"
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

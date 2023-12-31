import sys
import numpy as np
import os

hmm_clustering_path = os.path.join(os.getcwd(), "hmm_clustering")
sys.path.insert(0, hmm_clustering_path)
dataset_path = os.path.join(os.getcwd(), "dataset")
sys.path.insert(0, dataset_path)
from hmm_clustering.motifdata import Mtf
import hmm_clustering.greedycluster as greedy
import hmm_clustering.hmmmerge as hmm
import dataset.dbm as dbm
import dataset.dataset_yeast_mini as yst
import dataset.fungi_mini as fg
import dataset.fungi as fg2


# EVERYTHING IN THE ORDER OF A-C-G-T!
from sklearn.metrics import (
    fowlkes_mallows_score,
    adjusted_rand_score,
    normalized_mutual_info_score,
)


# def read_labels_from_file(file_path):
#     with open(file_path, "r") as f:
#         labels = [int(label) for label in f.read().split()]
#     return labels


""" 
Range: [0, 1]
# Interpretation: A value of 1 indicates perfect agreement and a value of 0
# indicates no agreement. 
"""


def evaluate_clustering(true_labels, predicted_labels):
    fmi = fowlkes_mallows_score(true_labels, predicted_labels)
    print("Fowlkes-Mallows Index (FMI):", fmi)


def evaluate_clustering_nmi(true_labels, predicted_labels):
    nmi = normalized_mutual_info_score(true_labels, predicted_labels)

    print("Normalized Mutual Information (NMI):", nmi)


def evaluate_clustering_ari(true_labels, predicted_labels):
    ari = adjusted_rand_score(true_labels, predicted_labels)
    print("Adjusted Rand Index (ARI):", ari)


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


# true_labels = read_labels_from_file("true_labels.txt")
# predicted_labels = read_labels_from_file("predicted_labels.txt")

# # Evaluate the clustering results
# evaluate_clustering(true_labels, predicted_labels)

# Dataset: fungi_mini
dataset2 = fg2.DNABindingMotifs()
# dataset2.mmm is a dictionary with format {"id": Bio.motifs.Motif object}

print(dataset2.mmm)
print("Finish printing dataset2 \n")

dataset2_cc = dataset2.cc
print(dataset2_cc)
print("Finish printing dataset2 cc \n")


# Convert dataset2_cc into a list of lists, where each inner list contains the keys of the corresponding dictionary
dataset2_cc_list = [[key for key in d] for d in dataset2_cc]
print(dataset2_cc_list)
print("Finish printing dataset2 cc list \n")

# Hierarchical clustering
# clusters = hierarchical_clustering.hierarchical_clustering(5, dm, idlist)
# # resulting cluters after using clustering method
# print(clusters)

num_clusters2 = len(dataset2_cc_list)  # Define the number of clusters
print(num_clusters2)
print("Finish printing num_clusters2 \n")

# distance_update can be one of ["complete","single","UPGMA","WPGMA"]
clusters2 = hmm.hmm_fin()
# clusters2 = kmeans_self_defined_dist.kmeans_motifs(dm2, idlist2, num_clusters2)
print(clusters2)
print("Finish printing result by HMM clustering \n")


# our_cluster2 = create_cluster_dicts(clusters2, dataset2)
# base_dir_our2 = "predicted_clusters2"
# generator_our2 = graph.MotifGraphGenerator(our_cluster2, base_dir_our2)
# cluster_paths_our2 = generator_our2.generate_cluster_graphs()
# generator_our2.generate_final_composite_graph(cluster_paths_our2)

print(dataset2_cc_list)
dataset2_dic = dict(sorted(generate_labels_from_clusters(dataset2_cc_list).items()))
print(dataset2_dic)
print("Finish printing true labels \n")

print(clusters2)
clusters2_dic = dict(sorted(generate_labels_from_clusters(clusters2).items()))
print(clusters2_dic)
print("Finish printing predicted labels \n")

true_labels2 = list(dataset2_dic.values())
predicted_labels2 = list(clusters2_dic.values())

score_fmi2 = evaluate_clustering(true_labels2, predicted_labels2)
print(score_fmi2)
score_nmi2 = evaluate_clustering_nmi(true_labels2, predicted_labels2)
print(score_nmi2)
score_ari2 = evaluate_clustering_ari(true_labels2, predicted_labels2)
print(score_ari2)

print(clusters2)
print("Finish printing result by HMM clustering \n")

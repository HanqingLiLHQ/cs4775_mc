import sys

sys.path.insert(1, "../CS4775_MC/dataset")
import dbm
import dataset_yeast_mini
import fungi_mini
from clustering import hierarchical_clustering
from clustering import kmeans, kmeans_self_defined_dist
from distance_matrix_calculation import pearson_distance
from distance_matrix_calculation import jensen_distance
from evaluate import ari_score
from evaluate import fmi_score
from evaluate import nmi_score

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


################################################################################

# Dataset: yeast_mini

dataset = dataset_yeast_mini.DNABindingMotifs()
# dataset.mmm is a dictionary with format {"id": Bio.motifs.Motif object}
# where each motif matrix can be retrieved by motif.pssm

print(dataset.mmm)
print("Finish printing dataset \n")

dm, idlist = pearson_distance.calculate_distance_matrix(dataset)
# dm, idlist = jensen_distance.calculate_distance_matrix(dataset)

dataset_cc = dataset.cc
print(dataset_cc)
print("Finish printing dataset cc \n")

# Convert dataset_cc into a list of lists, where each inner list contains the keys of the corresponding dictionary
dataset_cc_list = [[key for key in d] for d in dataset_cc]
print(dataset_cc_list)
print("Finish printing dataset cc list \n")

# K-means clustering
# One advantage of K-means is it can customize number of clusters needed
num_clusters = len(dataset_cc_list)  # Define the number of clusters
print(num_clusters)
print("Finish printing num_clusters \n")

clusters = kmeans_self_defined_dist.kmeans_motifs(dm, idlist, num_clusters)
print(clusters)
print("Finish printing result by K-means \n")

dataset_dic = dict(sorted(generate_labels_from_clusters(dataset_cc_list).items()))
print(dataset_dic)
print("Finish printing true labels \n")

clusters_dic = dict(sorted(generate_labels_from_clusters(clusters).items()))
print(clusters_dic)
print("Finish printing predicted labels \n")

true_labels = list(dataset_dic.values())
predicted_labels = list(clusters_dic.values())

score_ari = ari_score.evaluate_clustering(true_labels, predicted_labels)
score_fmi = fmi_score.evaluate_clustering(true_labels, predicted_labels)
score_nmi = nmi_score.evaluate_clustering(true_labels, predicted_labels)


################################################################################

# Dataset: fungi_mini

dataset2 = fungi_mini.DNABindingMotifs()
# dataset2.mmm is a dictionary with format {"id": Bio.motifs.Motif object}

print(dataset2.mmm)
print("Finish printing dataset2 \n")

dm2, idlist2 = pearson_distance.calculate_distance_matrix(dataset2)

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

clusters2 = kmeans_self_defined_dist.kmeans_motifs(dm2, idlist2, num_clusters2)
print(clusters2)
print("Finish printing result by K-means \n")

dataset2_dic = dict(sorted(generate_labels_from_clusters(dataset2_cc_list).items()))
print(dataset2_dic)
print("Finish printing true labels \n")

clusters2_dic = dict(sorted(generate_labels_from_clusters(clusters2).items()))
print(clusters2_dic)
print("Finish printing predicted labels \n")

true_labels2 = list(dataset2_dic.values())
predicted_labels2 = list(clusters2_dic.values())

score_ari2 = ari_score.evaluate_clustering(true_labels2, predicted_labels2)
score_fmi2 = fmi_score.evaluate_clustering(true_labels2, predicted_labels2)
score_nmi2 = nmi_score.evaluate_clustering(true_labels2, predicted_labels2)

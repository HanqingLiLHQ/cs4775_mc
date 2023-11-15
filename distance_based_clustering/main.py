import sys

sys.path.insert(1, "../CS4775_MC-1/dataset")
import dbm
import dataset_yeast_mini
from clustering import hierarchical_clustering
from clustering import kmeans, kmeans_self_defined_dist
from distance_matrix_calculation import pearson_distance
from distance_matrix_calculation import jensen_distance
from evaluate import ari_score
from evaluate import fmi_score
from evaluate import nmi_score


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


dataset_dic = dict(sorted(assign_unique_numbers_to_values(dataset.mmm).items()))
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

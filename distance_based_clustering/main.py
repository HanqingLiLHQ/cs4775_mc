import sys

sys.path.insert(1, "../CS4775_MC/dataset")
sys.path.insert(2, "../CS4775_MC/")
import dbm
import dataset_yeast_mini
import fungi_mini
import fungi
from column_distance import *
from clustering import hierarchical_clustering
from clustering import kmeans_self_defined_dist
from distance_matrix_calculation import pearson_distance
from distance_matrix_calculation import jensen_distance
from distance_matrix_calculation import dist_from_col
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

# Dataset: yeast_mini

dataset = dataset_yeast_mini.DNABindingMotifs()
# dataset.mmm is a dictionary with format {"id": Bio.motifs.Motif object}
# where each motif matrix can be retrieved by motif.pssm

print(dataset.mmm)
print("Finish printing dataset \n")

base_dir_motifs = "motifs_graphs"
motif_generator = motif_graph.MotifGraphGenerator(dataset.mmm, base_dir_motifs)
motif_generator.generate_motif_graphs()

dm, idlist = dist_from_col.calculate_distance_matrix(
    dataset, Euclidean_Distance, "expand", bg=[0.25, 0.25, 0.25, 0.25], average=np.mean
)
# dm, idlist = jensen_distance.calculate_distance_matrix(dataset)

dataset_cc = dataset.cc
print(dataset_cc)
print("Finish printing dataset cc \n")

# Convert dataset_cc into a list of lists, where each inner list contains the keys of the corresponding dictionary
dataset_cc_list = [[key for key in d] for d in dataset_cc]
print(dataset_cc_list)
print("Finish printing dataset cc list \n")

# base_dir = "nature_clusters"
# generator = graph.MotifGraphGenerator(dataset_cc, base_dir)
# cluster_paths = generator.generate_cluster_graphs()
# generator.generate_final_composite_graph(cluster_paths)

# For clusters of motifs
base_dir_clusters = "nature_clusters"
cluster_generator = cluster_graph.ClusterGraphGenerator(
    dataset_cc_list, base_dir_motifs, base_dir_clusters
)
cluster_paths = cluster_generator.generate_cluster_graphs()
final_graph.generate_final_composite_graph(base_dir_clusters, cluster_paths)

# K-means clustering
# One advantage of K-means is it can customize number of clusters needed
num_clusters = len(dataset_cc_list)  # Define the number of clusters
print(num_clusters)
print("Finish printing num_clusters \n")

# distance_update can be one of ["complete","single","UPGMA","WPGMA"]
clusters = hierarchical_clustering.hierarchical_clustering(
    num_clusters, dm, idlist, distance_update="complete"
)
print(clusters)
print("Finish printing result by Hierarchical Clustering \n")

# our_cluster = create_cluster_dicts(clusters, dataset)

# print(our_cluster)
# print("Finish printing our cluster \n")

base_dir_clusters_our = "predicted_clusters"
cluster_generator_our = cluster_graph.ClusterGraphGenerator(
    clusters, base_dir_motifs, base_dir_clusters_our
)
cluster_paths_our = cluster_generator_our.generate_cluster_graphs()

# base_dir_our = "predicted_clusters"
# generator_our = graph.MotifGraphGenerator(our_cluster, base_dir_our)
# cluster_paths_our = generator_our.generate_cluster_graphs()
# generator_our.generate_final_composite_graph(cluster_paths_our)


# clusters = kmeans_self_defined_dist.kmeans_motifs(dm, idlist, num_clusters)


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

base_dir_motifs2 = "motifs_graphs2"
motif_generator2 = motif_graph.MotifGraphGenerator(dataset2.mmm, base_dir_motifs2)
motif_generator2.generate_motif_graphs()

dm2, idlist2 = dist_from_col.calculate_distance_matrix(
    dataset2, Euclidean_Distance, "expand", bg=[0.25, 0.25, 0.25, 0.25], average=np.mean
)

dataset2_cc = dataset2.cc
print(dataset2_cc)
print("Finish printing dataset2 cc \n")

# base_dir2 = "nature_clusters2"
# generator2 = graph.MotifGraphGenerator(dataset2_cc, base_dir2)
# cluster_paths2 = generator2.generate_cluster_graphs()
# generator2.generate_final_composite_graph(cluster_paths2)


# Convert dataset2_cc into a list of lists, where each inner list contains the keys of the corresponding dictionary
dataset2_cc_list = [[key for key in d] for d in dataset2_cc]
print(dataset2_cc_list)
print("Finish printing dataset2 cc list \n")

# For clusters of motifs
base_dir_clusters2 = "nature_clusters2"
cluster_generator2 = cluster_graph.ClusterGraphGenerator(
    dataset2_cc_list, base_dir_motifs2, base_dir_clusters2
)
cluster_paths2 = cluster_generator2.generate_cluster_graphs()
final_graph.generate_final_composite_graph(base_dir_clusters2, cluster_paths2)


# Hierarchical clustering
# clusters = hierarchical_clustering.hierarchical_clustering(5, dm, idlist)
# # resulting cluters after using clustering method
# print(clusters)

num_clusters2 = len(dataset2_cc_list)  # Define the number of clusters
print(num_clusters2)
print("Finish printing num_clusters2 \n")

# distance_update can be one of ["complete","single","UPGMA","WPGMA"]
clusters2 = hierarchical_clustering.hierarchical_clustering(
    num_clusters2, dm2, idlist2, distance_update="complete"
)
# clusters2 = kmeans_self_defined_dist.kmeans_motifs(dm2, idlist2, num_clusters2)
print(clusters2)
print("Finish printing result by Hierarchical Clustering \n")

base_dir_clusters_our2 = "predicted_clusters2"
cluster_generator_our2 = cluster_graph.ClusterGraphGenerator(
    clusters2, base_dir_motifs2, base_dir_clusters_our2
)
cluster_paths_our2 = cluster_generator_our2.generate_cluster_graphs()
final_graph.generate_final_composite_graph(base_dir_clusters_our2, cluster_paths_our2)


# our_cluster2 = create_cluster_dicts(clusters2, dataset2)
# base_dir_our2 = "predicted_clusters2"
# generator_our2 = graph.MotifGraphGenerator(our_cluster2, base_dir_our2)
# cluster_paths_our2 = generator_our2.generate_cluster_graphs()
# generator_our2.generate_final_composite_graph(cluster_paths_our2)


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

################################################################################

# Dataset: fungi

dataset3 = fungi.DNABindingMotifs()
# dataset3.mmm is a dictionary with format {"id": Bio.motifs.Motif object}

print(dataset3.mmm)
print("Finish printing dataset3 \n")

base_dir_motifs3 = "motifs_graphs3"
motif_generator3 = motif_graph.MotifGraphGenerator(dataset3.mmm, base_dir_motifs3)
motif_generator3.generate_motif_graphs()


dm3, idlist3 = dist_from_col.calculate_distance_matrix(
    dataset3, Euclidean_Distance, "expand", bg=[0.27, 0.23, 0.23, 0.27], average=np.mean
)

dataset3_cc = dataset3.cc
print(dataset3_cc)
print("Finish printing dataset3 cc \n")

# base_dir3 = "nature_clusters3"
# generator3 = graph.MotifGraphGenerator(dataset3_cc, base_dir3)
# cluster_paths3 = generator3.generate_cluster_graphs()
# generator3.generate_final_composite_graph(cluster_paths3)

# Convert dataset3_cc into a list of lists, where each inner list contains the keys of the corresponding dictionary
dataset3_cc_list = [[key for key in d] for d in dataset3_cc]
dataset3_cc_list.sort(key=len, reverse=True)
# dataset3_cc_list.reverse()
dataset3_cc_list = sort_sublists(dataset3_cc_list)
print(dataset3_cc_list)
print("Finish printing dataset3 cc list \n")

# For clusters of motifs
base_dir_clusters3 = "nature_clusters3"
cluster_generator3 = cluster_graph.ClusterGraphGenerator(
    dataset3_cc_list, base_dir_motifs3, base_dir_clusters3
)
cluster_paths3 = cluster_generator3.generate_cluster_graphs()
final_graph.generate_final_composite_graph(base_dir_clusters3, cluster_paths3)


num_clusters3 = len(dataset3_cc_list)  # Define the number of clusters
print(num_clusters3)
print("Finish printing num_clusters3 \n")

# distance_update can be one of ["complete","single","UPGMA","WPGMA"]
cluster3 = hierarchical_clustering.hierarchical_clustering(
    num_clusters3, dm3, idlist3, distance_update="complete"
)
cluster3 = sort_sublists(cluster3)
print(cluster3)
print("Finish printing result by hierarchical clustering \n")

cluster3_match = most_similar_lists(dataset3_cc_list, cluster3)
print(cluster3_match)
print("Finish printing result by matching \n")

base_dir_clusters_our3 = "predicted_clusters3"
cluster_generator_our3 = cluster_graph.ClusterGraphGenerator(
    cluster3_match, base_dir_motifs3, base_dir_clusters_our3
)
cluster_paths_our3 = cluster_generator_our3.generate_cluster_graphs()
final_graph.generate_final_composite_graph(base_dir_clusters_our3, cluster_paths_our3)

# cluster2 = kmeans_self_defined_dist.kmeans_motifs(dm2, idlist2, num_clusters2)
# print(cluster2)
# print("Finish printing result by K means \n")

# our_cluster3 = create_cluster_dicts(cluster3, dataset3)

# base_dir_our3 = "predicted_clusters3"
# generator_our3 = graph.MotifGraphGenerator(our_cluster3, base_dir_our3)
# cluster_paths_our3 = generator_our3.generate_cluster_graphs()
# generator_our3.generate_final_composite_graph(cluster_paths_our3)


dataset3_dic = dict(sorted(generate_labels_from_clusters(dataset3_cc_list).items()))
print(dataset3_dic)
print("Finish printing true labels \n")

clusters3_dic = dict(sorted(generate_labels_from_clusters(cluster3).items()))
print(clusters3_dic)
print("Finish printing predicted labels \n")

true_labels3 = list(dataset3_dic.values())
predicted_labels3 = list(clusters3_dic.values())

score_ari3 = ari_score.evaluate_clustering(true_labels3, predicted_labels3)
score_fmi3 = fmi_score.evaluate_clustering(true_labels3, predicted_labels3)
score_nmi3 = nmi_score.evaluate_clustering(true_labels3, predicted_labels3)

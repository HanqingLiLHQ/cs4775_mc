import sys
import numpy as np

sys.path.insert(1, "../CS4775_MC")  # Add the parent directory to sys.path

from hmm_clustering.motifdata import Mtf 
import hmm_clustering.greedycluster as greedy
import hmm_clustering.hmmmerge as hmm
import dataset.dbm as dbm
import dataset.dataset_yeast_mini as yst
import dataset.fungi_mini as fg
import dataset.fungi as fg2
import main as scorer


#EVERYTHING IN THE ORDER OF A-C-G-T!

from sklearn.metrics import (
    fowlkes_mallows_score,
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
'''
  Initializes data and returns motifs stored in Mtf object 

  returns:
  mtf_dict: dictionaries of all the Motifs in Mtf format, with key=id
    and value=Mtf object
'''
def init_mtf():
    mtf_dict = {}
    curr_id = 0
    db_instance = fg2.DNABindingMotifs().mmm
    # print(db_instance)
    for i in db_instance:
        temp_mtf = Mtf(mtf=db_instance[i], id=curr_id, label = i)
        # pwm is a dictionary; j.pwm.get('A') gives back a very long tuple.
        my_id = temp_mtf.get_id()
        # print(temp_mtf.get_pwm())
        # print(temp_mtf.get_cols())
        mtf_dict[my_id] = temp_mtf
        curr_id += 1
    gddct = greedy.greedyclus(0.3, mtf_dict)
    # print(gddct.keys())
    #v1 = [[0.1,0.2,0.3,0.4],[0.1,0.2,0.3,0.4],[0.12,0.18,0.25,0.45],[0.1,0.2,0.3,0.4]]
    #v2 = [[0.1,0.2,0.3,0.4],[0.1,0.2,0.3,0.4],[0.26,0.24,0.25,0.25],[0.12,0.18,0.25,0.45],[0.1,0.2,0.3,0.4]]
    #v1a = np.array(v1)
    #v2a = np.array(v2)
    #print(hmm.hmm_scoring_align(v1a, v2a))


    for i in range(23):
        if len(gddct) == 2:
            break
        gddct = hmm.hmm_merge(gddct)
        print(gddct.keys())


    motif_list = []
    mtf_list = [value[2] for value in gddct.values()]
    for mtfclus in mtf_list:
        motif_clus = [mtf.get_label() for mtf in mtfclus]
        motif_list = motif_list + [motif_clus]
    #print(motif_list)
    return motif_list


def main():

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
    clusters2 = init_mtf()
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


main()
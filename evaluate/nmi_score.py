from sklearn.metrics import (
    normalized_mutual_info_score,
)


# def read_labels_from_file(file_path):
#     with open(file_path, "r") as f:
#         labels = [int(label) for label in f.read().split()]
#     return labels


"""
Range: [0, 1]
Interpretation: A value of 1 implies perfect agreement between the
two clusterings, while a value of 0 implies no agreement.
"""


def evaluate_clustering(true_labels, predicted_labels):
    nmi = normalized_mutual_info_score(true_labels, predicted_labels)

    print("Normalized Mutual Information (NMI):", nmi)


# true_labels = read_labels_from_file("true_labels.txt")
# predicted_labels = read_labels_from_file("predicted_labels.txt")

# # Evaluate the clustering results
# evaluate_clustering(true_labels, predicted_labels)

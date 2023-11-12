from sklearn.metrics import (
    adjusted_rand_score,
)


def read_labels_from_file(file_path):
    with open(file_path, "r") as f:
        labels = [int(label) for label in f.read().split()]
    return labels


"""
Range: [-1, 1]
Interpretation: A value of 1 implies that the two sets of clusters are
identical (i.e., the predicted labels perfectly match the true labels).
A value close to 0 implies that the clusters are assigned randomly.
A value less than 0 indicates less than random agreement.
"""


def evaluate_clustering(true_labels, predicted_labels):
    ari = adjusted_rand_score(true_labels, predicted_labels)
    print("Adjusted Rand Index (ARI):", ari)


true_labels = read_labels_from_file("true_labels.txt")
predicted_labels = read_labels_from_file("predicted_labels.txt")

# Evaluate the clustering results
evaluate_clustering(true_labels, predicted_labels)

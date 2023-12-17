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


# true_labels = read_labels_from_file("true_labels.txt")
# predicted_labels = read_labels_from_file("predicted_labels.txt")

# # Evaluate the clustering results
# evaluate_clustering(true_labels, predicted_labels)

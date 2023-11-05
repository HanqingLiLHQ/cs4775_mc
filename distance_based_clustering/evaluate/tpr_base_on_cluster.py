def calculate_tpr(cluster, ground_truth):
    """
    Calculate the True Positive Rate for a single cluster against the ground truth.

    param cluster: a set of motif instances that are clustered together
    param ground_truth: a set of motif instances that are known to be true(correct labeling)
    return: TPR as a float
    """
    tp = len(cluster.intersection(ground_truth))  # True positives
    fn = len(ground_truth) - tp  # False negatives

    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0
    return tpr


# Assuming each cluster and the ground_truth are sets of motif instance identifiers
clusters = [
    {"motif1", "motif2", "motif3"},
    {"motif4", "motif5"},
    {"motif6", "motif7"},
]

# This would be your known ground truth data
ground_truth = {"motif1", "motif2", "motif3", "motif6", "motif7"}

# Calculate TPR for each cluster
tprs = [calculate_tpr(cluster, ground_truth) for cluster in clusters]

# Now tprs holds the TPR for each cluster
print(tprs)

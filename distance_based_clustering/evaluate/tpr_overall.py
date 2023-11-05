def calculate_tpr(clusters, ground_truth):
    """
    Calculate the overall True Positive Rate for all clusters against the ground truth.

    param clusters: a list of sets, each set contains motif instances that are clustered together
    param ground_truth: a set of motif instances that are known to be true
    return: overall TPR as a float
    """
    total_tp = 0  # Total true positives
    total_fn = 0  # Total false negatives

    # Go through each cluster and count TPs and FNs
    for cluster in clusters:
        # True positives in the current cluster
        tp = len(cluster.intersection(ground_truth))

        # False negatives in the current cluster
        fn = len(ground_truth.difference(cluster))

        total_tp += tp
        total_fn += fn

    # Adjust for double-counted FNs across clusters
    total_fn -= (len(clusters) - 1) * len(ground_truth)

    # Calculate overall TPR
    overall_tpr = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    return overall_tpr


# Example usage:
clusters = [
    {"motif1", "motif2", "motif3"},
    {"motif2", "motif4"},
    # ... more clusters
]

ground_truth = {"motif1", "motif2", "motif3", "motif4", "motif5"}

overall_tpr = calculate_tpr(clusters, ground_truth)
print(f"Overall TPR: {overall_tpr}")

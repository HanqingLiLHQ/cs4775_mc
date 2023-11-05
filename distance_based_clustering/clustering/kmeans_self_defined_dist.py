import random
import numpy as np

def select_initial_centroids(num_clusters, dm):
    # Randomly choose 'num_clusters' different indices to serve as initial centroids
    indices = list(range(len(dm)))
    random.shuffle(indices)
    centroids_indices = indices[:num_clusters]
    return centroids_indices

def assign_clusters(centroids_indices, dm):
    clusters = [[] for _ in range(len(centroids_indices))]
    for i, distances_to_centroids in enumerate(dm):
        # Find the closest centroid for each motif
        closest_centroid = np.argmin([dm[i][c] for c in centroids_indices])
        clusters[closest_centroid].append(i)
    return clusters

def find_new_centroids(clusters, dm):
    new_centroids_indices = []
    for cluster in clusters:
        # Find the motif in the cluster that has the minimum average distance to all other motifs in the cluster
        if cluster:
            min_avg_distance = float('inf')
            centroid_index = None
            for motif_index in cluster:
                avg_distance = np.mean([dm[motif_index][other_index] for other_index in cluster])
                if avg_distance < min_avg_distance:
                    min_avg_distance = avg_distance
                    centroid_index = motif_index
            new_centroids_indices.append(centroid_index)
    return new_centroids_indices

def kmeans_motifs(dm, id_list, num_clusters, max_iter=100):
    # Initialize centroids
    centroids_indices = select_initial_centroids(num_clusters, dm)
    
    for iteration in range(max_iter):
        # Assign clusters based on the current centroids
        clusters = assign_clusters(centroids_indices, dm)
        
        # Find new centroids within the current clusters
        new_centroids_indices = find_new_centroids(clusters, dm)
        
        # Check for convergence (if centroids do not change)
        if set(new_centroids_indices) == set(centroids_indices):
            break
        
        centroids_indices = new_centroids_indices

    # Transform the indices to the actual motif ids
    clusters_with_ids = [[id_list[index] for index in cluster] for cluster in clusters]
    return clusters_with_ids
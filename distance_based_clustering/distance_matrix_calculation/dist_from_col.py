from Bio import motifs
import sys

sys.path.insert(1, "../CS4775_MC")
import motif_distance
from  column_distance import *
import numpy as np

def calculate_distance_matrix(dataset):
    # Create a PWM for each motif and normalize it with pseudocounts
    ppm_map = {
        motif_id: motif.counts.normalize(pseudocounts=0.5) 
        for motif_id, motif in dataset.mmm.items()
    }

    # Extract motif IDs for indexing
    motif_ids = list(ppm_map.keys())

    # Initialize the distance matrix with zeros
    distance_matrix = np.zeros((len(motif_ids), len(motif_ids)))

    # Calculate the Jensen-Shannon Divergence for each pair of motifs
    for i, motif_id_1 in enumerate(motif_ids):
        for j, motif_id_2 in enumerate(motif_ids):
            if i < j:  # Since the matrix is symmetric, we only need to calculate once
                ppm1 = np.array(list(ppm_map[motif_id_1].values())).T
                ppm2 = np.array(list(ppm_map[motif_id_2].values())).T
                # Kullback_Leibler_Distance
                # Jensen_Shannon_Distance
                # Euclidean_Distance
                distance = motif_distance.distance(Kullback_Leibler_Distance, ppm1, ppm2,"expand")
                distance_matrix[i][j] = distance_matrix[j][i] = distance

    return distance_matrix.tolist(), motif_ids
from Bio import motifs
from scipy.spatial.distance import jensenshannon
import numpy as np

def calculate_distance_matrix(dataset):
    # Create a PWM for each motif and normalize it with pseudocounts
    pwm_map = {
        motif_id: motif.counts.normalize(pseudocounts=0.5) 
        for motif_id, motif in dataset.mmm.items()
    }

    # Extract motif IDs for indexing
    motif_ids = list(pwm_map.keys())

    # Initialize the distance matrix with zeros
    distance_matrix = np.zeros((len(motif_ids), len(motif_ids)))

    # Calculate the Jensen-Shannon Divergence for each pair of motifs
    for i, motif_id_1 in enumerate(motif_ids):
        for j, motif_id_2 in enumerate(motif_ids):
            if i < j:  # Since the matrix is symmetric, we only need to calculate once
                pwm1 = np.array(list(pwm_map[motif_id_1].values()))
                pwm2 = np.array(list(pwm_map[motif_id_2].values()))
                # Calculate JSD for each column and take the average
                distance = np.mean([jensenshannon(col1, col2) for col1, col2 in zip(pwm1.T, pwm2.T)])
                distance_matrix[i][j] = distance_matrix[j][i] = distance

    return distance_matrix.tolist(), motif_ids
from Bio import motifs
import sys

sys.path.insert(1, "../CS4775_MC")
import motif_distance
from  column_distance import *
import numpy as np
from scipy import stats

def calculate_distance_matrix(dataset, col_dist, align_method, bg = [0.25,0.25,0.25,0.25], average = np.mean):
    """
    Return the distance matrix of the whole dataset, calculated via column-wise distance measurement [col_dist] with [align_method]. 
    If the alignment method is chosen to be "expand", the algorithm will calculate a threshold to see whether the minimum value
    really indicates a meaningful motif alignment. If not, the algorithm will return the mean distance of all possible alignments (as all are possible). 
    The column comparison scores are averaged via [average] function. 
    Preconditions:
        col_dist: returns a numerical distance value based on two input 4D vectors
            The current implementation supports: 
                Kullback_Leibler_Distance, Jensen_Shannon_Distance, Euclidean_Distance, Pearson_CC_Distance,
            with Pearson_CC_Distance only available when bg is set not to be [0.25,0.25,0.25,0.25]. 
        bg: background frequency, represented by 4-D vectors that should sum to 1 (only expand_compare will use it)
        alignment_method: ["expand","overlap"]
        average: some average method that return a numerical value with a input list
    """
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
                # The first attribute could be:
                # Kullback_Leibler_Distance
                # Jensen_Shannon_Distance
                # Euclidean_Distance
                # average method could be:

                distance = motif_distance.distance(col_dist, ppm1, ppm2,align_method, bg = bg, average = average)
                distance_matrix[i][j] = distance_matrix[j][i] = distance

    return distance_matrix.tolist(), motif_ids
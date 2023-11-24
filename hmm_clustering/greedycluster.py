# Sort sequences
# Initialize first sequence as a cluster
# Starting from the second sequence, each sequence is compared to all current clusters.
# If the sequence meet a pre-defined threshold for a current cluster, it joins that cluster.
# If the sequence does not meet a pre-defined threshold, it starts a new representative cluster.

import numpy as np
from pyjaspar import jaspardb
import sys

sys.path.insert(1, "../CS4775_MC")  # Add the parent directory to sys.path

import column_distance as coldis
import motif_distance as mtfdis
from hmm_clustering.motifdata import Mtf  

#EVERYTHING IN THE ORDER OF A-C-G-T!

'''
  Computes the distance between two mtf objects

  returns:
  dist: float representing the calculated distance
'''
def eulc_mtf_distcalc(mtf1: Mtf, mtf2: Mtf):
  ppm1 = mtf1.get_cols()
  ppm2 = mtf2.get_cols()
  return mtfdis.distance(coldis.Euclidean_Distance, ppm1, ppm2, "expand")

def greedyclus(mtf_dict):
    results = []

    keys = list(mtf_dict.keys())

    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            id1, id2 = keys[i], keys[j]
            mtf1, mtf2 = mtf_dict[id1], mtf_dict[id2]

            distance = eulc_mtf_distcalc(mtf1, mtf2)

            result_str = f"id1={id1}, id2={id2}, dist={distance}"
            print(result_str)
            results.append(result_str)

    return "\n".join(results)
    






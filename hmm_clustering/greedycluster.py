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
  Computes the distance between two mtf objects using the euclidian method.

  Input:
  mtf1, mtf2: Mtf objects
  returns:
  dist: float representing the calculated distance
'''
def eulc_mtf_distcalc(mtf1: Mtf, mtf2: Mtf):
  ppm1 = mtf1.get_cols()
  ppm2 = mtf2.get_cols()
  return mtfdis.distance(coldis.Euclidean_Distance, ppm1, ppm2, "expand")

'''
TOCHANGE
  Computes the distance between all pairs of matrices and print them out.

  Input:
  mtf_dict: dictionaries of all the Motifs in Mtf format, with key=id
    and value=Mtf object
  returns:
  TOCHANGE
'''
def greedyclustemp(mtf_dict):
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

'''
TOCHANGE
  Computes the distance between all pairs of matrices and print them out.

  Input:
  mtf_dict: dictionaries of all the Motifs in Mtf format, with key=id
    and value=Mtf object
  returns:
  TOCHANGE
'''
def greedyclus(thresh, mtf_dict):
  assert len(mtf_dict) != 0
  threshold = thresh
  clusters = {} #finalized clusters, dict of int-(cols, [ids], [Mtfs])
  clusters[0] = (mtf_dict[0].get_cols(), [0], [mtf_dict[0]])
  curr_mtx_id = 1 # id of the next matrix to compare
  curr_clus_id = 1 #id of the next cluster
  keys = list(mtf_dict.keys())
  # loop through the matrices
  while curr_mtx_id < len(keys):
    clus = -1
    min_dist = -1
    curr_clus = 0
    # test distance on the new matrix candidate to each existing clusters
    for (cols, ids, mtfs) in clusters.values():
      temp_dist = mtfdis.distance(coldis.Euclidean_Distance, cols, mtf_dict[curr_mtx_id].get_cols(), "expand")
      # attempts to find minimum distance smaller than threshold value
      if temp_dist <= threshold and (clus == -1 or temp_dist <= min_dist):
        clus = curr_clus
        min_dist = temp_dist
      curr_clus += 1

    # no distance under threshold was found-create new cluster
    if clus == -1:
       clusters[curr_clus_id] = (mtf_dict[curr_mtx_id].get_cols(), [curr_mtx_id], [mtf_dict[curr_mtx_id]])
       curr_clus_id += 1
    # adds new Mtf object to the identified cluster and recalculate columns
    else:
       target_cluster = clusters[clus]
       modified_ids = target_cluster[1] + [curr_mtx_id]
       modified_mtfs = target_cluster[2] + [mtf_dict[curr_mtx_id]]
       modified_cols = calculate_mean_vectors([mymtf.get_cols() for mymtf in modified_mtfs])
       clusters[clus] = (modified_cols, modified_ids, modified_mtfs)
    curr_mtx_id += 1
  print(clusters)
  return clusters

"""
Helper function for updating the pwm of multiple Mtf objects.

Input:
cols_list: list of 4*n matrix
Output:
newcols: 4*n matrix of updated means
"""
def calculate_mean_vectors(cols_list):
    newcols_list = []
    # Find the maximum length of vectors in the list
    max_length = max(len(cols) for cols in cols_list)
    # Pad vectors with zeros to make them of equal length
    padded_colist = []
    for cols in cols_list: 
      num_rows_to_append = max_length - len(cols)
      pm = np.pad(cols, ((0, num_rows_to_append), (0, 0)), mode='constant', constant_values=0)
      padded_colist.append(pm)
      padded_cols = np.array(padded_colist)
    # print(padded_cols)
    for i in range(max_length):
      cols_array = []
      for padcols in padded_cols:
          if any(padcols[i] != 0):
            cols_array.append(padcols[i])

      newcols_list.append(np.mean(cols_array, axis=0))
    newcols = np.array(newcols_list)
    # print("newcols = "+str(newcols))
    return newcols


       









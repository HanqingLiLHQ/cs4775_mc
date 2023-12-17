import numpy as np
from pyjaspar import jaspardb
import sys

sys.path.insert(1, "../CS4775_MC")  # Add the parent directory to sys.path

import column_distance as coldis
import motif_distance as mtfdis
from hmm_clustering.motifdata import Mtf  
import hmm_clustering.greedycluster as greedy
import dataset.dbm as dbm
import dataset.dataset_yeast_mini as yst
import dataset.fungi_mini as fg
import dataset.fungi as fg2

#EVERYTHING IN THE ORDER OF A-C-G-T!


"""
Helper function for aligning the new input motif according to given offset.

Input:
offset: Offset Position? Yes-int(>=0)
ppm1: first mtf cols for input
ppm2: second mtf cols for input
Output:
fixedmtfs: list of Mtf objects with fixated cols
"""
def offsetchop(offset, ppm1, ppm2):
   assert len(ppm1) >= len(ppm2) + offset
   return ppm1[offset:(len(ppm2) + offset)], ppm2

"""
  Resort two ppms to put the longer one in the front.
Input:
  ppm1: first ppm
  ppm2: second ppm
Output:
  longerppm, shorterppm: longer ppm in the first output position, shorter one 
in the second output position.
"""
def arrange_ppm(ppm1, ppm2):
   if len(ppm1)>=len(ppm2):
      return ppm1, ppm2
   else:
      return ppm2, ppm1


"""
Stuff [0,0,0,0]s to the shorter ppm until they are of the same length.
Input:
  ppm1: first ppm
  ppm2: second ppm
Output:
  ppm1: first ppm
  ppm2: second ppm
"""
def paddle_ppms(ppm1, ppm2):
    max_length = max(len(ppm1), len(ppm2))
    # Pad vectors with zeros to make them of equal length
    padded_colist = []
    for cols in [ppm1, ppm2]: 
      num_rows_to_append = max_length - len(cols)
      pm = np.pad(cols, ((0, num_rows_to_append), (0, 0)), mode='constant', constant_values=0)
      padded_colist.append(pm)
      padded_cols = np.array(padded_colist)
    return padded_colist[0], padded_colist[1]

"""
  Returns the best alignment of the representative cluster PPMs, 
list of on/off targets during matching,  as well as the calculated 
scores based on the constructed HMM models. The HMM models of each
cluster is similar to a profile-HMM model, but given the preconditions of
ungapped motifs, the Insertion and Deletion states are replaced by an 
"offtarget" state with emmission probabilities of each of the four nucleotides
as 0.25. 
Input:
  ppm1 : representative motif of the first cluster, np.array object
  ppm2 : representative motif of the second cluster, np.array object
Precondition:
  ppm1 and ppm2 has the same length
Output:
  offset: integer representing the offset from mtf1 to mtf2
  l: list of most likely hidden states at each position, on/off target
  s: score from the list of hidden states upon best alignment
"""
def hmm_score(ppm1, ppm2):
  assert len(ppm1) == len(ppm2)
  # Ontarget/Offtarget
  n_states = 2
  n_obs = len(ppm1)
  # equal prob of on/off target
  trans_probs = {0:{0:0.5, 1:0.5}, 1:{0:0.5,1:0.5}}
  states = list(trans_probs.keys())
  # Initialize matrices
  # Row 0 is ontarget, Row 1 is offtarget
  vb_m = np.zeros((n_states, n_obs))
  btr_m = np.zeros((n_states, n_obs), dtype=int)

  # Initialize the first column of the vbm, using log likelihood to calculate emission probabilities
  vb_m[0][0] = coldis.Euclidean_Distance(ppm1[0],ppm2[0])
  vb_m[1][0] = coldis.Euclidean_Distance([0.25,0.25,0.25,0.25],ppm2[0])

  # Recurse through the rest of the vbm
  for t in range(1, n_obs):
      for curr_state in range(n_states):
          # Find max aikV(t-1)
          min_score = float("inf")
          best_prev_state = None

          for prev_state in range(n_states):
              score = vb_m[prev_state][t - 1]

              if score < min_score:
                  min_score = score
                  best_prev_state = prev_state
          # Multiply by ek(xt)
          if curr_state == 0:
             vb_m[curr_state][t] = min_score + coldis.Euclidean_Distance(ppm1[t],ppm2[t])
          elif curr_state == 1:
             vb_m[curr_state][t] = min_score + coldis.Euclidean_Distance([0.25,0.25,0.25,0.25],ppm2[t])
          # unmatched
          else:
             vb_m[curr_state][t] = min_score

          btr_m[curr_state][t] = best_prev_state

              # Backtracking
    
  # final state:[obs[-1]]
  l = [0] * n_obs
  best_final_state = np.argmin(vb_m[:, -1])

  for t in range(n_obs, 0, -1):
      t -= 1  
      l[t] = best_final_state
      best_final_state = btr_m[best_final_state][t]
  # penalty for offmatching is 0.5
  p = vb_m[best_final_state][-1] + sum(l) * 0.5

  return l, p
  

"""
  Conducting position-wise HMM-scoring between two ppms, giving out the minimum 
  score, the best offset and the according on/off target sequence. 
Input:
  ppm1 : representative motif of the first cluster, np.array object
  ppm2 : representative motif of the second cluster, np.array object
Output:
  minscore: minimum score after position-wise comparison
  bestoffset: best offset for alignment (shorter matrix to longer matrix)
  l: list of best HMM states at best alignment
"""
def hmm_scoring_align(ppm1, ppm2):
   (longppm, shortppm) = arrange_ppm(ppm1, ppm2)
   offset = len(longppm) - len(shortppm)
   l = None
   minscore = float("inf")
   bestoffset = -1
   for i in range(offset + 1):
      (tempppm1, tempppm2) = offsetchop(i, longppm, shortppm)
      (myalignment, myscore) = hmm_score(tempppm1, tempppm2)
      if myscore < minscore:
         minscore = myscore
         bestoffset = i
         l = myalignment
   return minscore, bestoffset, l


"""
  Merge two given clusters.
Input:
  clusters: naive clusters from greedyclustering, 
    dictionary of int-(cols, [ids], [Mtfs])
  offset: offset for merging alignments.
  i1: the first key of the cluster for merging.
  i2: the second key of the cluster for merging.
Output:
  merged_clusters: merged clusters, dictionary of 
  int-(cols, [ids], [Mtfs])
"""
def merge(clusters, offset, i1, i2):
    # Check if i1 and i2 are valid keys in the clusters dictionary
    if i1 not in clusters or i2 not in clusters:
        print("Invalid keys for merging.")
        return clusters

    # Extract information for the two clusters to be merged
    cols1, ids1, mtfs1 = clusters[i1]
    cols2, ids2, mtfs2 = clusters[i2]

    # Merge information (you can customize this based on your requirements)
    merged_ids = ids1 + ids2
    merged_mtfs = offsetfix(offset, mtfs1, mtfs2)
    merged_cols = greedy.calculate_mean_vectors([mymtf.get_cols() for mymtf in merged_mtfs])

    # Create the merged cluster
    merged_cluster = (merged_cols, merged_ids, merged_mtfs)

    # Remove the old clusters from the dictionary
    del clusters[i1]
    del clusters[i2]

    # Add the merged cluster to the dictionary with a new key (you may want to customize this key)
    new_key = max(clusters.keys()) + 1
    clusters[new_key] = merged_cluster

    return clusters

"""
  Conduct step-wise merge based on _scoring_align.
Input:
  clusters: naive clusters from greedyclustering, 
    dictionary of int-(cols, [ids], [Mtfs])
Output:
  clusters: finalized clusters for scoring
"""
def hmm_merge(clusters):
   best_minscore = float("inf")
   best_offset = None
   best_l = []
   key_list = list(clusters.keys())
   ppm_list = [value[0] for value in clusters.values()]
   id1 = -1
   id2 = -1
   for i in range(len(ppm_list)):
      for j in range(i + 1, len(ppm_list)):
         minscore, offset, l = hmm_scoring_align(ppm_list[i], ppm_list[j])
         if minscore < best_minscore:
            best_minscore = minscore
            best_l = l
            best_offset = offset
            id1 = key_list[i]
            id2 = key_list[j]
   print("merging clusters with id " + str(id1)+" and "+str(id2)+" with score "+str(best_minscore)+" and HMM states "+str(best_l))
   return merge(clusters, best_offset, id1, id2)

"""
Helper function for aligning the new input motif according to given offset.

Input:
offset: Offset Position? Yes-int, No-None
mtfs1: list of Mtf objects
mtfs2: list of Mtf objects
Output:
fixedmtfs: list of Mtf objects with fixated cols
"""
def offsetfix(offset, mtfs1, mtfs2):
   if offset == None or offset == 0:
      return mtfs1 + mtfs2
   # [0, len(motif1)) and [offset, len(motif2) + offset)]
   elif offset > 1:
      for tempmtf in mtfs2:
         prealigned_cols = tempmtf.get_cols()
         zeros_array = np.zeros_like(prealigned_cols[:1])
         aligned_cols = np.vstack([zeros_array] * offset * -1 + [prealigned_cols])
         tempmtf.set_cols(aligned_cols)
      return mtfs1 + mtfs2
   else:
      for tempmtf in mtfs1:
         prealigned_cols = tempmtf.get_cols()
         zeros_array = np.zeros_like(prealigned_cols[:1])
         aligned_cols = np.vstack([zeros_array] * offset * -1 + [prealigned_cols])
         tempmtf.set_cols(aligned_cols)
      return mtfs1 + mtfs2
   

'''
  Initializes data and returns motifs stored in Mtf object 

  returns:
  mtf_dict: dictionaries of all the Motifs in Mtf format, with key=id
    and value=Mtf object
'''
def hmm_fin():
    mtf_dict = {}
    curr_id = 0
    db_instance = fg2.DNABindingMotifs().mmm
    # print(db_instance)
    for i in db_instance:
        temp_mtf = Mtf(mtf=db_instance[i], id=curr_id, label = i)
        # pwm is a dictionary; j.pwm.get('A') gives back a very long tuple.
        my_id = temp_mtf.get_id()
        # print(temp_mtf.get_pwm())
        # print(temp_mtf.get_cols())
        mtf_dict[my_id] = temp_mtf
        curr_id += 1
    gddct = greedy.greedyclus(0, mtf_dict)
    # print(gddct.keys())
    #v1 = [[0.1,0.2,0.3,0.4],[0.1,0.2,0.3,0.4],[0.12,0.18,0.25,0.45],[0.1,0.2,0.3,0.4]]
    #v2 = [[0.1,0.2,0.3,0.4],[0.1,0.2,0.3,0.4],[0.26,0.24,0.25,0.25],[0.12,0.18,0.25,0.45],[0.1,0.2,0.3,0.4]]
    #v1a = np.array(v1)
    #v2a = np.array(v2)
    #print(hmm.hmm_scoring_align(v1a, v2a))


    for i in range(113):
        if len(gddct) == 2:
            break
        gddct = hmm_merge(gddct)
        print(gddct.keys())


    motif_list = []
    mtf_list = [value[2] for value in gddct.values()]
    for mtfclus in mtf_list:
        motif_clus = [mtf.get_label() for mtf in mtfclus]
        motif_list = motif_list + [motif_clus]
    #print(motif_list)
    return motif_list

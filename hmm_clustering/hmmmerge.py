import numpy as np
from pyjaspar import jaspardb
import sys

sys.path.insert(1, "../CS4775_MC")  # Add the parent directory to sys.path

import column_distance as coldis
import motif_distance as mtfdis
from hmm_clustering.motifdata import Mtf  
import hmm_clustering.greedycluster as greedy

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
   alignment = None
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

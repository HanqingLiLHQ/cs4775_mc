import sys
sys.path.insert(1,'../../dataset')
import dbm
from Bio import motifs



"""
Return a motif-to-motif distance matrix as a 2D list, with a 1D list composed of motif ids 
to indicate the motif's order in the distance matrix
"""
def calculate_distance_matrix():
  # gather the dataset, currently using the DNABindingMotifs dataset
  dataset = dbm.DNABindingMotifs()
  # gather the cluster

  mmm = dataset.mmm
  motif_ids = list(mmm.keys())

  # pseudocount all the motifs before calculating the PSSM (position specific scoring matrix)
  for motif_id in motif_ids:
    mmm[motif_id].pseudocounts = motifs.jaspar.calculate_pseudocounts(mmm[motif_id])
     
  # generate a 2D matrix that denotes a distance matrix. The correspondance of the matrix 
  # index and the motif_ids is denoted with a list. 
  mmd_matrix = []
  for i in range(len(motif_ids)):
    motif1 = motif_ids[i]
    motif1_distances = []
    for j in range(len(motif_ids)):
      motif2 = motif_ids[j]
      motif1_distances.append((mmm[motif1].pssm).dist_pearson(mmm[motif2].pssm)[0])
    mmd_matrix.append(motif1_distances)

  return mmd_matrix, motif_ids

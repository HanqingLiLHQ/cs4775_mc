from Bio import motifs


"""
Return a motif-to-motif distance matrix as a 2D list, with a 1D list composed of motif ids 
to indicate the motif's order in the distance matrix
"""


def calculate_distance_matrix(dataset):
    mmm = dataset.mmm
    motif_ids = list(mmm.keys())

    # pseudocount all the motifs before calculating the PSSM (position specific scoring matrix)
    # Computes the root square of the total number of sequences multiplied by the background nucleotide.
    # the background nucleotide frequency is uniform if not specified.
    for motif_id in motif_ids:
        mmm[motif_id].pseudocounts = motifs.jaspar.calculate_pseudocounts(
            mmm[motif_id])

    # generate a 2D matrix that denotes a distance matrix. The correspondance of the matrix
    # index and the motif_ids is denoted with a list.
    # pssm matrix is calculated log_2(P/P_background) for each position after the pseudocount is added.
    mmd_matrix = []
    for i in range(len(motif_ids)):
        motif1 = motif_ids[i]
        motif1_distances = []
        for j in range(len(motif_ids)):
            motif2 = motif_ids[j]
            motif1_distances.append(
                (mmm[motif1].pssm).dist_pearson(mmm[motif2].pssm)[0])
        mmd_matrix.append(motif1_distances)

    return mmd_matrix, motif_ids


# The calculation of pearson distance is as follows:
# calculate the pearson score for the ungapped spanning of one motif through the other,
# find the largest possible PCC score, and return 1 - PCC. The offset is defined by where the
# other motif should be aligned.

# Specically, the pearson score is calculated only for the overlapping part, but the normalization
# factor includes the total length of the alignment. Therefore, a larger PCC would prefer more overlapping.
# Calculate the total sum of values in each matrix as sx and sy.
# calculate the total sum of squared values in each aligned matrix as sxx, syy
# calculate the sum of pairwise products in the two matrix as sxy
# normalize them by the total amount of grids (total alignment length * possible bases)
# PCC = (sxy - sx * sy) / root((sxx - sx * sx) (syy - sy * sy))
# which is just: PCC = Cov(mx, my)/ root(Var(x) * Var(y))

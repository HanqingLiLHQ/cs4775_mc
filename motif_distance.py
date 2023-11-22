# This module uses column comparison methods from column_distance to compare between motifs


from column_distance import *
import numpy as np
from Bio import motifs
from Bio.motifs import Motif
from pyjaspar import jaspardb
import matplotlib.pyplot as  plt

def preprocess(motif):
    """
    Return pseudocounted ppm of a motif object
    """
    return np.array(list(motif.counts.normalize(pseudocounts=0.1).values())).T


def naive_compare(distance_method,ppm1, ppm2):
    """
    This version compare the columns one by one, both starting from the first. The column distance similarities are averaged
    to calculate the motif similarity. 
    Precondition:
        distance_method: returns a numerical distance value based on two input 4D vectors
        ppm1, ppm2: n * 4 matrices with the same dimensionality
    """
    assert ppm1.shape == ppm2.shape
    col_dists = [distance_method(col1, col2) for col1, col2 in zip(ppm1, ppm2)]
    distance = np.mean(col_dists)
    return distance
    

def add_column(ppm, start, end, col=[0.25,0.25,0.25,0.25]):
    """
    helper function to add flanking position at backgrouund frequency
    Precondition:
        ppm: n * 4 numpy matrix, each row represents a position
        start, end: non-negative integers indicating how many new positions to add to the ppm at the [start/end]
        col: background frequency, represented by 4-D vectors that should sum to 1 
    """
    start_m = np.array(start * [col])
    end_m = np.array(end * [col])
    if start == end == 0:
        return ppm
    if start == 0:
        return np.vstack([ppm,end_m])
    elif end == 0:
        return np.vstack([start_m, ppm])
    return np.vstack([start_m, ppm,end_m])


def expand_compare(col_dist,ppm1,ppm2,bg=[0.25,0.25,0.25,0.25]):
    """
    Slide one motif through the other to check all possible ungapped alignments, the unmatched positions are supplied with
    background probability.
    Preconditions:
        col_dist: returns a numerical distance value based on two input 4D vectors
        ppm1/ppm2: n * 4 numpy matrix, each row represents a position
        bg: background frequency, represented by 4-D vectors that should sum to 1 
    """
    off_dist = {}
    # swap the two matrices to ensure that ppm1 is longer than ppm2
    if len(ppm2) > len(ppm1):
        ppm2, ppm1 = ppm1, ppm2
    len1 = len(ppm1)
    len2 = len(ppm2)
    assert len1 >= len2
    # iterate through all possible offsets
    for offset in range(-len(ppm2)+1,len(ppm1)):
        if offset < 0: 
            exp_ppm2 = add_column(ppm2, 0, len1 - offset - len2, col = bg)
            exp_ppm1 = add_column(ppm1, -offset, 0, col = bg)
        elif offset >= 0 and offset <= len1 - len2:
            exp_ppm2 = add_column(ppm2, offset, len1 - len2 - offset, col = bg)
            exp_ppm1 = ppm1
        elif offset > len1 - len2:
            exp_ppm2 = add_column(ppm2, offset, 0, col = bg)
            exp_ppm1 = add_column(ppm1, 0, offset + len2 - len1, col = bg)
        off_dist[offset] = naive_compare(col_dist, exp_ppm1, exp_ppm2)
    return off_dist
            
def cut_compare(col_dist, ppm1, ppm2):
    """
    Slide one motif through the other to check all possible ungapped alignments. Only consider the overlapping region
    Preconditions:
        col_dist: returns a numerical distance value based on two input 4D vectors
        ppm1/ppm2: n * 4 numpy matrix, each row represents a position
    """    
    off_dist = {}
    # swap the two matrices to ensure that ppm1 is longer than ppm2
    if len(ppm2) > len(ppm1):
        ppm2, ppm1 = ppm1, ppm2
    len1 = len(ppm1)
    len2 = len(ppm2)
    assert len1 >= len2
    for offset in range(-len(ppm2)+1,len(ppm1)):
        if offset < 0: 
            cut_ppm2 = ppm2[-offset:]
            cut_ppm1 = ppm1[:(len2 + offset)]
        elif offset >= 0 and offset <= len1 - len2:
            cut_ppm2 = ppm2
            cut_ppm1 = ppm1[offset:len2 + offset]
        elif offset > len1 - len2:
            cut_ppm1 = ppm1[offset:]
            cut_ppm2 = ppm2[:len1 - offset]
        off_dist[offset] = naive_compare(col_dist, cut_ppm1, cut_ppm2)
    return off_dist



def compare_align(col_dist,ppm1,ppm2, bg = [0.25,0.25,0.25,0.25]):
    """
    Helper function to compare the expand/cut motif comparison methods given two ppm and the distance method
    Preconditions:
        col_dist: returns a numerical distance value based on two input 4D vectors
        ppm1/ppm2: n * 4 numpy matrix, each row represents a position
        bg: background frequency, represented by 4-D vectors that should sum to 1 (only expand_compare will use it)
    """
    plot, ax = plt.subplots()
    exp_results = expand_compare(col_dist, ppm1, ppm2,bg = bg)
    cut_results = cut_compare(col_dist, ppm1, ppm2)
    plt.plot(exp_results.keys(), exp_results.values(),label = "Expanded Motif Comparison")
    plt.plot(cut_results.keys(), cut_results.values(), label = "Overlap Motif Comparison")
    plt.legend()
    plt.title("Comparison of Two Motif Alignment (sliding) Methods")
    ax.set_xticks(list(exp_results.keys()))
    ax.set_ylim(bottom = 0, top = None)
    ax.set_xlabel("offset (positions)")
    ax.set_ylabel("motif distance")
    return


def compare_distance(ppm1,ppm2, bg = [0.25,0.25,0.25,0.25]):
    """
    Helper function to compare different column-wise distance measurements in terms of sliding alignment of expanded motifs
    Preconditions:
        ppm1/ppm2: n * 4 numpy matrix, each row represents a position
    """
    plot, ax = plt.subplots()
    l2_results = expand_compare(Euclidean_Distance, ppm1, ppm2, bg = bg)
    js_results = expand_compare(Jensen_Shannon_Distance,ppm1,ppm2,bg = bg)
    kl_results = expand_compare(Kullback_Leibler_Distance,ppm1,ppm2, bg = bg)
    l2_values = np.array(list(l2_results.values())) / sum(l2_results.values())
    js_values = np.array(list(js_results.values())) / sum(js_results.values())
    kl_values = np.array(list(kl_results.values())) / sum(kl_results.values())
    plt.plot(l2_results.keys(),l2_values,label = "Euclidean_Distance")
    plt.plot(js_results.keys(),js_values,label = "Jensen_Shannon_Distance")
    plt.plot(kl_results.keys(), kl_values,label = "Kullback_Leibler_Distance")
    plt.title("Comparison of Different Column-wise Distance Measurements with Expanded Motif Alignment")
    plt.legend()
    ax.set_xticks(list(l2_results.keys()))
    ax.set_ylim(bottom = 0, top = None)
    ax.set_xlabel("offset (positions)")
    ax.set_ylabel("motif distance")
    return


def distance(col_dist, ppm1, ppm2, align_method, bg = [0.25,0.25,0.25,0.25]):
    """
    Return the shortest alignment distance between ppm1, ppm2, calculated via column-wise distance measurement [col_dist] with [align_method]. 
    Preconditions:
        col_dist: returns a numerical distance value based on two input 4D vectors
        ppm1/ppm2: n * 4 numpy matrix, each row represents a position
        bg: background frequency, represented by 4-D vectors that should sum to 1 (only expand_compare will use it)
        alignment_method: ["expand","overlap"]
    """
    if align_method == "expand":
        return min(list(expand_compare(col_dist, ppm1,ppm2, bg = bg).values()))
    elif align_method == "overlap":
        return min(list(cut_compare(col_dist, ppm1, ppm2).values()))
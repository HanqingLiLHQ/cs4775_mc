# Sort sequences
# Initialize first sequence as a cluster
# Starting from the second sequence, each sequence is compared to all current clusters.
# If the sequence meet a pre-defined threshold for a current cluster, it joins that cluster.
# If the sequence does not meet a pre-defined threshold, it starts a new representative cluster.

import numpy as np
from pyjaspar import jaspardb
from Bio import motifs
from Bio.motifs import Motif
import motifdata as mtf
import sys

sys.path.insert(1, "../CS4775_MC")  # Add the parent directory to sys.path

import column_distance as coldis
from hmm_clustering.motifdata import Mtf  


def mtf_distcalc(mtflst):
  pass





# Sort sequences
# Initialize first sequence as a cluster
# Starting from the second sequence, each sequence is compared to all current clusters.
# If the sequence meet a pre-defined threshold for a current cluster, it joins that cluster.
# If the sequence does not meet a pre-defined threshold, it starts a new representative cluster.

from dataset.dataset_yeast_mini import DNABindingMotifs
from pyjaspar import jaspardb
from Bio import motifs
from Bio.motifs import Motif
from hmm_clustering.motifdata import Mtf
import random


'''
  Initializes data and returns the pwm stored in a list. 
  pwms are power density matrices that stored as dictionaries with 
    keys 'A', 'T', 'C', 'G' and values as tuples.

  returns:
  pwm_list: list of pwms from the source data
'''


def init_mtf():
    mtf_list = []
    curr_id = 0
    db_instance = DNABindingMotifs()
    for i in db_instance.dbs:
        for j in i:
            temp_mtf = Mtf(mtf=j, id=curr_id)
            # pwm is a dictionary; j.pwm.get('A') gives back a very long tuple.
            print(temp_mtf.id)
            mtf_list.append(temp_mtf)
            curr_id += 1
    return mtf_list


def main():
    init_mtf()


main()

import sys
import numpy as np

sys.path.insert(1, "../CS4775_MC")  # Add the parent directory to sys.path

from hmm_clustering.motifdata import Mtf 
import hmm_clustering.greedycluster as greedy
import hmm_clustering.hmmmerge as hmm
import dataset.dbm as dbm
import dataset.dataset_yeast_mini as yst
import dataset.fungi_mini as fg
import dataset.fungi as fg2


#EVERYTHING IN THE ORDER OF A-C-G-T!
'''
  Initializes data and returns motifs stored in Mtf object 

  returns:
  mtf_dict: dictionaries of all the Motifs in Mtf format, with key=id
    and value=Mtf object
'''
def init_mtf():
    mtf_dict = {}
    curr_id = 0
    db_instance = fg2.DNABindingMotifs()
    for i in db_instance.dbs:
        for j in i:
            temp_mtf = Mtf(mtf=j, id=curr_id)
            # pwm is a dictionary; j.pwm.get('A') gives back a very long tuple.
            my_id = temp_mtf.get_id()
            # print(temp_mtf.get_pwm())
            # print(temp_mtf.get_cols())
            mtf_dict[my_id] = temp_mtf
            curr_id += 1
    # greedy.greedyclus(0.2, mtf_dict)
    v1 = [[0.1,0.2,0.3,0.4],[0.1,0.2,0.3,0.4],[0.12,0.18,0.25,0.45],[0.1,0.2,0.3,0.4]]
    v2 = [[0.1,0.2,0.3,0.4],[0.1,0.2,0.3,0.4],[0.26,0.24,0.25,0.25],[0.12,0.18,0.25,0.45],[0.1,0.2,0.3,0.4]]
    v1a = np.array(v1)
    v2a = np.array(v2)
    print(hmm.hmm_scoring_align(v1a, v2a))
    return mtf_dict


def main():
    init_mtf()


main()
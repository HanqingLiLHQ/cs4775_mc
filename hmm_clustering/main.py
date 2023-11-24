import sys

sys.path.insert(1, "../CS4775_MC")  # Add the parent directory to sys.path

from hmm_clustering.motifdata import Mtf 
import hmm_clustering.greedycluster as greedy
import dataset.dbm as dbm
import dataset.dataset_yeast_mini as yst
import dataset.fungi_mini as fg


'''
  Initializes data and returns the pwm stored in a list. 
  pwms are power density matrices that stored as dictionaries with 
    keys 'A', 'T', 'C', 'G' and values as tuples.

  returns:
  pwm_list: list of pwms from the source data
'''
def init_mtf():
    mtf_dict = {}
    curr_id = 0
    db_instance = fg.DNABindingMotifs()
    for i in db_instance.dbs:
        for j in i:
            temp_mtf = Mtf(mtf=j, id=curr_id)
            # pwm is a dictionary; j.pwm.get('A') gives back a very long tuple.
            my_id = temp_mtf.get_id()
            mtf_dict[my_id] = temp_mtf
            curr_id += 1
    print(mtf_dict)
    return mtf_dict


def main():
    init_mtf()


main()
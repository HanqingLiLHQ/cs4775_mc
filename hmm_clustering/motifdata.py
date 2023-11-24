from dataset.dataset_yeast_mini import DNABindingMotifs
from pyjaspar import jaspardb
from Bio import motifs
from Bio.motifs import Motif
import random
"""
  Custom class that contains a Motif data to be analyzed

  mtf: Motif object
  id: int that stores a unique id for each value in a dataset
"""


class Mtf():
    def __init__(self, mtf, id):
        self.mtf = mtf
        self.id = id


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
  mtf = None
  id = None
  def __init__(self, mtf, id):
      self.mtf = mtf
      self.id = id
    
  """
  Returns power weight matrix of a motif object

  returns:
  pwm: pwm is a dictionary; pwm.get('A') gives back a very long tuple.
  """
  def get_pwm(self):
      return self.mtf.pwm
  
  """
  Returns id of a Mtf object

  returns:
  id: an int representing a unique Mtf object.
  """
  def get_id(self):
      return self.id
  

  """
  Prints details of a Mtf object.
  """
  def __str__(self):
      return "TODO"


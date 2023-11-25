
from pyjaspar import jaspardb
from Bio import motifs
from Bio.motifs import Motif
import numpy as np

#EVERYTHING IN THE ORDER OF A-C-G-T!

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
  pwms are power weight matrices that stored as dictionaries with 
    keys 'A', 'C', 'G', 'T' and values as tuples.

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
  Returns all the columns of the counts of a Mtf object in a list

  returns:
  cols_normalized: n*4 matrix, with each of the the vectors having probabilities sum 
    up to 1, that represents the probability of A,C,G,T respectively.
  """
  def get_cols(self):
    a = self.mtf.counts.get('A')
    t = self.mtf.counts.get('T')
    c = self.mtf.counts.get('C')
    g = self.mtf.counts.get('G')

    # Assuming all a, t, c, g have the same length
    cols = np.array([a, c, g, t]).T  # Transpose to get columns

    cols_normalized = np.apply_along_axis(lambda x: x / np.sum(x), axis=1, arr=cols)

    return cols_normalized
  
    """
  Returns all the columns of the pwm of a Mtf object in a list

  returns:
  cols: n*4 matrix, with each of the the vectors having probabilities sum 
    up to 1, that represents the probability of A,C,G,T respectively.
  """
  def get_pwm_cols(self):
    a = self.mtf.pwm.get('A')
    t = self.mtf.pwm.get('T')
    c = self.mtf.pwm.get('C')
    g = self.mtf.pwm.get('G')

    # Assuming all a, t, c, g have the same length
    cols = np.array([a, c, g, t]).T  # Transpose to get columns

    return cols

  """
  Returns the ppm as a n*4 list of a Mtf object

  returns:
  ppm: 4*n matrix representing the pwm of a Mtf object
  """
  def get_ppm(self):
    a = self.mtf.pwm.get('A')
    c = self.mtf.pwm.get('C')
    g = self.mtf.pwm.get('G')
    t = self.mtf.pwm.get('T')

    ppm = np.vstack([a, c, g, t])
    
    return ppm
  

  """
  Prints details of a Mtf object.
  """
  def __str__(self):
      return "TODO"


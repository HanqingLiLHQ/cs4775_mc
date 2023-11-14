# This file contains functions that output a distance given two columns of a
# positional count matrix.

# Most of the methods come from this paper: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-2-r24
import numpy as np


def Euclidean_Distance(col1, col2):
    """
      Evaluate the Euclidean Distance between two columns of the probability weight matrix (column sum = 1)
      col1, col2: a 4-dimension vector, the probabilities sum up to 1
    """
    col1 = np.array(col1)
    col2 = np.array(col2)
    # Precondition check
    assert col1.sum() == 1
    assert col2.sum() == 1
    return np.sqrt(((col1 - col2) * (col1 - col2)).sum())


def Pearson_CC_Distance(col1, col2):
    """
      Evaluate the Pearson correlation coefficient between two columns of the probability weight matrix (column sum = 1)
      Return 1 - PCC as a distance measurement. 
      The 0.25 in the implementation referes to the average value of the four entries in a vector
      col1, col2: a 4-dimension vector, the probabilities sum up to 1
    """
    col1 = np.array(col1)
    col2 = np.array(col2)
    # Precondition check
    assert col1.sum() == 1
    assert col2.sum() == 1
    numerator = np.sum((col1 - 0.25) * (col2 - 0.25))
    denominator = np.sqrt(np.sum((col1 - 0.25) ** 2)
                          * np.sum((col2 - 0.25) ** 2))
    return 1 - numerator / denominator


def Kullback_Leibler_Distance(col1, col2):
    """
        Evaluate the Kullback_Leibler_Divergence between two columns of the probability weight matrix (column sum = 1)
        As the KLD is not a symmetric distance measurement, a modified version that calculates the divergences for both directions
        and averages them is implemented here. 
        Precondition:
            col1, col2: a 4-dimension vector with non-zero probability, the probabilities sum up to 1
    """
    col1 = np.array(col1)
    col2 = np.array(col2)
    # Precondition check
    assert col1.sum() == 1
    assert col2.sum() == 1
    assert np.array(col1 == 0).sum() == 0
    assert np.array(col2 == 0).sum() == 0
    return 0.5 * (np.sum(col1 * np.log(col1 / col2)) + np.sum(col2 * np.log(col2 / col1)))
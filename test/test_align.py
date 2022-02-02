# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    test = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat",-10,-1)
    test.align("MQR", "MYQR")                     
                   

    ## Calculate this matrices manually
    ## Looking only at the first two rows
    manual_alignmatrix = np.matrix([[0,-(np.inf),-(np.inf),-(np.inf),-(np.inf)], [-(np.inf), 5,-12,-12,-14]])
    manual_gapA = np.matrix([[-10,-(np.inf),-(np.inf),-(np.inf),-(np.inf)], [-(11), -22,-23,-24,-25]])  
    manual_gapB = np.matrix([[-10,-11,-12,-13,-14], [-(np.inf), -22,-6,-7,-8]])  

    assert(np.all(test._align_matrix[:2] == manual_alignmatrix))
    assert(np.all(test._gapA_matrix[:2] == manual_gapA))
    assert(np.all(test._gapB_matrix[:2] == manual_gapB))



def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    test = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat",-10,-1)
    print(test.align("MQR", "MYQR"))




test_nw_backtrace()
test_nw_alignment()
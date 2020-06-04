from scipy.stats import entropy
from RCMAP.Alignment import Alignments
from RCMAP.Evaluation_seq import *

def frequence_aa(aa,pos):
    if  alignments.count_aa_ref()[pos-1][aa] == 0 :
        return 0
    return alignments.count_aa_ref()[pos-1][aa]/len(alignments.seqrefs)

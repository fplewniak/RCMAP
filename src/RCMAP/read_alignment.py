from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition

seq = list(SeqIO.parse('ArsM_aln.faa', 'fasta'))
print (seq)
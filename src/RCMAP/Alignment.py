from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from RCMAP.Classification_AA import AAcategories


class Alignments :

    def __init__(self, file ="ArsM_aln.faa", seqs_to_evaluate =["WP_045226361.1", "Q969Z2"]):
        self.file = file
        self.seqs_to_evaluate = seqs_to_evaluate
        alignment = AlignIO.read(file, "fasta")
        self.seqeval = MultipleSeqAlignment([s for s in alignment if s.id in seqs_to_evaluate])
        self.seqrefs = MultipleSeqAlignment([s for s in alignment if s.id not in seqs_to_evaluate])

    def read_alignments(self):
        return self.seqrefs, self.seqeval

    def get_cat_at_pos(self, pos):
        AA_at_pos = set()
        for k in range(len(self.seqrefs)) :
            AA_at_pos.add(self.seqrefs[k][pos])
        return (AAcategories().find_category(AA_at_pos))

    #def get_cat_in_range(seqref,debut=None,fin=None):
        #return [get_cat_at_pos(seqref,pos) for pos in range(debut,fin)]

object = Alignments()
#print(object.read_alignments(file,seqs_to_evaluate))
print(object.get_cat_at_pos(3))
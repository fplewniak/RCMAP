from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from RCMAP.Classification_AA import AAcategories


class Alignments :

    def __init__(self, file, seqs_to_evaluate):
        self.file = file
        self.seqs_to_evaluate = seqs_to_evaluate
        alignment = AlignIO.read(file, "fasta")
        self.seqeval = MultipleSeqAlignment([s for s in alignment if s.id in seqs_to_evaluate])
        self.seqrefs = MultipleSeqAlignment([s for s in alignment if s.id not in seqs_to_evaluate])

    def get_alignments(self):
        return self.seqrefs, self.seqeval

    def get_cat_at_pos(self, pos):
        """
        :param pos: position of the amino acid
        :return:
        """
        AA_at_pos = set()
        for k in range(len(self.seqrefs)) :
            AA_at_pos.add(self.seqrefs[k][pos])
        return (AAcategories().find_category(AA_at_pos))


    def get_cat_in_range(self):
        """
        :param pos1: start position #from 1 until end
        :param pos2: end position #from 1 until end
        :return: list of the categories of amino acid at every position between pos1 and pos2
        """
        cat_in_range = []
        for k in range(len(self.get_positions_list(positions))):
            if self.get_positions_list(positions)[k][0] == None :
                self.get_positions_list(positions)[k][0] = 0
            if self.get_positions_list(positions)[k][1] == None :
                self.get_positions_list(positions)[k][1] = len(self.seqrefs)
            l =[]
            for i in range(self.get_positions_list(positions)[k][0],self.get_positions_list(positions)[k][1]):
                l.append(self.get_cat_at_pos(i))
            cat_in_range.append(l)
        return cat_in_range

#object = Alignments("ArsM_aln.faa",["WP_045226361.1", "Q969Z2"])
#print(object.get_alignments(file,seqs_to_evaluate))
#print(object.get_cat_at_pos(2))
#print(object.get_cat_in_range(None,None))
positions =["0:3",":3","2:",":"]
print(object. get_positions_list(positions))
print(object.get_cat_in_range())




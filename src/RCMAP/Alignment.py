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


#    def get_cat_in_range(self,pos1=None,pos2=None):
#        """
#        :param pos1: start position #from 1 until end
#        :param pos2: end position #from 1 until end
#        :return: list of the categories of amino acid at every position between pos1 and pos2
#        """
#        if pos1 == None :
#            pos1=1
#        if pos2 == None :
#            pos2=len(self.seqrefs[0])
#        if pos1 > len(self.seqrefs[0]) or pos2 > len(self.seqrefs[0]):
#            return "Error"
#        cat_in_range = []
#        for pos in range(pos1-1, pos2):
#            cat_in_range.append(self.get_cat_at_pos(pos))
#        return cat_in_range

    def get_cat_in_range_bis(self, position):
        cat_in_range_bis = []
        if len(position) == 1 :
            cat_in_range_bis.append(self.get_cat_at_pos(position[0]))
        else:
            if position[0] == None:
                position[0]=1
            if position[1] == None :
                position[1]=len(self.seqrefs[0])
            if position[0] > len(self.seqrefs[0]) or position[1] > len(self.seqrefs[0]):
                return "Error"
            for pos in range(position[0]-1,position[1]):
                cat_in_range_bis.append(self.get_cat_at_pos(pos))
        return cat_in_range_bis

    def get_cat_in_range(self,interval):
        get_cat_in_range =[]
        if len(interval)==1:
            get_cat_in_range.append(self.get_cat_at_pos(interval[0]))
        else :
            for i in range(len(interval)):
                get_cat_in_range.append(self.get_cat_in_range_bis(interval[i]))
        return get_cat_in_range


object = Alignments("ArsM_aln.faa",["WP_045226361.1", "Q969Z2"])
#print(object.get_alignments(file,seqs_to_evaluate))
#print(object.get_cat_at_pos(2))
#print(object.get_cat_in_range(None,None))
print(object.get_cat_in_range([0,3]))
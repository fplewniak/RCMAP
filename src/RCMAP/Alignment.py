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
        :param pos: position of the amino acid in seqrefs
        :return:
        """
        AA_at_pos = set()
        for k in range(len(self.seqrefs)) :
            AA_at_pos.add(self.seqrefs[k][pos - 1])
        return (AAcategories().find_category(AA_at_pos))

    def get_aa_at_pos_in_seqeval(self, pos, name_seqeval):
        """
        :param pos: position of the amino acid in seqeval
        :return:
        """
        AA_at_pos_in_seqeval = set()
        for s in self.seqeval :
            if s.id == name_seqeval:
                AA_at_pos_in_seqeval = set(s[pos-1])
        return AA_at_pos_in_seqeval

    def get_aa_in_range_in_seqeval(self,pos1=None,pos2=None,name_seqeval):
        """
        :param pos1: start position #from 1 until end
        :param name_seqeval : name of the evaluated sequence
        :param pos2: end position #from 1 until end
        :return: list of amino acid at every position between pos1 and pos2 in seqeval
        """
        for s in range(self.seqeval):
            if s.id == name_seqeval:
                seqeval = s.seq
        if pos1 == None :
            pos1=1
        if pos2 == None :
            pos2=len(seqeval)
        #if pos1 > len(self.seqrefs[0]) or pos2 > len(self.seqrefs[0]):
        #    return "Error"
        aa_in_range = []
        for pos in range(pos1, pos2 + 1):
            aa_in_range.append(self.get_aa_at_pos_in_seqeval(pos,name_seqeval))
        return aa_in_range

    def get_cat_in_range(self,pos1=None,pos2=None):
        """
        :param pos1: start position #from 1 until end
        :param pos2: end position #from 1 until end
        :return: list of the categories of amino acid at every position between pos1 and pos2
        """
        if pos1 == None :
            pos1=1
        if pos2 == None :
            pos2=len(self.seqrefs[0])
        #if pos1 > len(self.seqrefs[0]) or pos2 > len(self.seqrefs[0]):
        #    return "Error"
        cat_in_range = []
        for pos in range(pos1, pos2 + 1):
            cat_in_range.append(self.get_cat_at_pos(pos))
        return cat_in_range

    def get_category_list(self,positions_list):
        """
        :param positions_list: list of the intervals of positions
        :return: list of the categories associated so the intervals of positions
        """
        list_of_categories =[]
        for r in range(len(positions_list)):
            if len(positions_list[r]) == 1:
                list_of_categories.append(self.get_cat_at_pos(positions_list[r][0]))
            else :
                list_of_categories.append(self.get_cat_in_range(positions_list[r][0],positions_list[r][1]))
        return list_of_categories


object = Alignments("ArsM_aln.faa",["WP_045226361.1", "Q969Z2"])
#print(object.get_alignments(file,seqs_to_evaluate))
#print(object.get_cat_at_pos(2))
#print(object.get_cat_in_range(None,None))
#positions =["0:3",":3","2:",":"]
#print(object. get_positions_list(positions))
#print(object.get_cat_in_range())




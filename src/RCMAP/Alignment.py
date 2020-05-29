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
        self.aa_ref_counts = self.count_aa_ref()
        self.list_of_categories = self.determine_ref_categories()

    def count_aa_ref(self):
        """
        :return: the count of amino acids at every position in all reference sequences
        """
        self.aa_ref_counts = [{"A":0,"R":0,"N":0,"D":0,"B":0,"C":0,"E":0,"Q":0,"Z":0,"G":0,"H":0,
                        "I":0,"L":0,"K":0,"M":0,"F":0,"P":0,"S":0,"T":0,"W":0,"Y":0,"V":0,"-":0} for sub in range(len(self.seqrefs[0]))]
        for s in self.seqrefs :
            for pos in range(len(self.seqrefs[0])):
                self.aa_ref_counts[pos][self.get_aa_at_pos(pos+1)]
        return self.aa_ref_counts

    def determine_ref_categories(self):
        """
        :return: the list of categories of amino acids at every position in seqrefs
        """
        self.list_of_categories = [set() for sub in range(len(self.seqrefs[0]))]
        for pos in range(len(self.aa_ref_counts)):
            self.list_of_categories[pos] = AAcategories().find_category({AA for AA in self.aa_ref_counts[pos] if self.aa_ref_counts[pos][AA] >0})

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
        return AAcategories().find_category(AA_at_pos)

    def get_aa_at_pos_in_seqeval(self, name_seqeval,pos):
        """
        :param pos: position of the amino acid in seqeval
        :return:
        """
        AA_at_pos_in_seqeval = set()
        for s in self.seqeval :
            if s.id == name_seqeval:
                AA_at_pos_in_seqeval = set(s[pos-1])
        return AA_at_pos_in_seqeval

    def get_aa_in_range_in_seqeval(self,name_seqeval,pos1=None,pos2=None):
        """
        :param pos1: start position #from 1 until end
        :param name_seqeval : name of the evaluated sequence
        :param pos2: end position #from 1 until end
        :return: list of amino acid at every position between pos1 and pos2 in seqeval
        """
        for s in self.seqeval:
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
            aa_in_range.append(self.get_aa_at_pos_in_seqeval(name_seqeval,pos))
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
                list_of_categories.append([self.get_cat_at_pos(positions_list[r][0])])
            else :
                list_of_categories.append(self.get_cat_in_range(positions_list[r][0],positions_list[r][1]))
        return list_of_categories

    def get_aa_list_in_seqeval(self,name_seqeval,positions_list) :
        """
        :param positions_list: list of the intervals of positions
        :return: list of the categories associated so the intervals of positions
        """
        aa_list_in_seqeval =[]
        for r in range(len(positions_list)):
            if len(positions_list[r]) == 1:
                aa_list_in_seqeval.append([self.get_aa_at_pos_in_seqeval(name_seqeval,positions_list[r][0])])
            else :
                aa_list_in_seqeval.append(self.get_aa_in_range_in_seqeval(name_seqeval,positions_list[r][0],positions_list[r][1]))
        return aa_list_in_seqeval






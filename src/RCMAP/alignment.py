from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from RCMAP.classification_aa import AAcategories
from scipy.stats import entropy


class Alignments:
    """
    Alignment is used to open the file and process the alignment of sequences
    """

    def __init__(self, file, seqs_to_evaluate):
        self.file = file
        self.seqs_to_evaluate = seqs_to_evaluate
        try:
            self.alignment = AlignIO.read(file, "fasta")
        except FileNotFoundError:
            print('File {fichier} not found'.format(fichier=file))
            exit(1)
        for seq in seqs_to_evaluate:
            k = 0
            for i in self.alignment:
                if i.id == seq:
                    k += 1
            if k == 0:
                print('Sequence to evaluate : {seq} not found'.format(seq=seq))
                exit(1)
        self.seqeval = MultipleSeqAlignment(
            [seq for seq in self.alignment if seq.id in seqs_to_evaluate])
        self.seqrefs = MultipleSeqAlignment(
            [s for s in self.alignment if s.id not in seqs_to_evaluate])
        self.aa_ref_counts = self.count_aa_ref()
        self.list_of_categories, self.list_of_cat_sets, self.list_of_aa_ref, self.count_aa_pos, \
        self.set_of_aa_ref = self.determine_ref_categories()

    def count_aa_ref(self):
        """
        :return: the count of amino acids at every position in all reference sequences
        """
        self.aa_ref_counts = [
            {"A": 0, "R": 0, "N": 0, "D": 0, "B": 0, "C": 0, "E": 0, "Q": 0, "Z": 0, "G": 0, "H": 0,
             "I": 0, "L": 0, "K": 0, "M": 0, "F": 0, "P": 0, "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0,
             "-": 0} for sub in range(len(self.seqrefs[0]))]
        for seq in self.seqrefs:
            for pos in range(len(self.seqrefs[0])):
                self.aa_ref_counts[pos][self.get_aa_at_pos(pos + 1, seq.id)] += 1
        return self.aa_ref_counts

    def determine_ref_categories(self):
        """
        :return: the list of categories of amino acids at every position in seqrefs
        """
        self.list_of_categories = [set() for sub in range(len(self.seqrefs[0]))]
        self.list_of_cat_sets = [set() for sub in range(len(self.seqrefs[0]))]
        self.list_of_aa_ref = []
        self.count_aa_pos = []
        self.set_of_aa_ref = []
        for pos in range(len(self.count_aa_ref())):
            self.list_of_categories[pos], self.list_of_cat_sets[pos] = \
                AAcategories().find_category(
                    {aa for aa in self.aa_ref_counts[pos] if self.aa_ref_counts[pos][aa] > 0})
            self.list_of_aa_ref.append([aa for aa in self.aa_ref_counts[pos] if
                                        self.aa_ref_counts[pos][aa] > 0])
            self.set_of_aa_ref.append({aa for aa in self.aa_ref_counts[pos] if
                                       self.aa_ref_counts[pos][aa] > 0})
            self.count_aa_pos.append(
                [self.aa_ref_counts[pos][aa] for aa in self.aa_ref_counts[pos] if
                 self.aa_ref_counts[pos][aa] > 0])
        return self.list_of_categories, self.list_of_cat_sets, self.list_of_aa_ref, \
               self.count_aa_pos, self.set_of_aa_ref

    def get_alignments(self):
        """
        :return: alignment of the reference sequences, alignment of the evaluated sequences
        """
        return self.seqrefs, self.seqeval

    def get_cat_set_at_pos(self, pos):
        """
        :param pos: position of the amino acid in seqrefs
        :return: the name of the category of amino acids observed in seqrefs
        """
        return self.list_of_cat_sets[pos - 1]

    def get_cat_at_pos(self, pos, strict):
        """
        :param strict: if True, the category is only the amino acids observed
        :param pos: position of the amino acid in seqrefs
        :return: the category (set) of amino acids observed in seqrefs
        """
        if strict:
            return self.set_of_aa_ref[pos - 1]
        return self.list_of_categories[pos - 1]

    def get_aa_observed_at_pos(self, pos):
        """
        :param pos:  position of the amino acids in seqrefs
        :return: all the amino acids observed in seqrefs at this position and the count of this
        amino acids
        """
        list_aa, list_count = self.list_of_aa_ref[pos - 1], self.count_aa_pos[pos - 1]
        aa_observed = dict()
        for i in range(len(list_aa)):
            aa_observed[list_aa[i]] = list_count[i]
        return dict(sorted(aa_observed.items(), key=lambda item: item[1], reverse=True))

    def get_aa_at_pos(self, pos, name_seq):
        """
        :param pos: position of the amino acid in seqref or seqeval
        :param name_seq: name of the sequence in seqref or seqeval
        :return: the amino acid at this position in the sequence
        """
        aa_at_pos = set()
        try:
            for seq in self.alignment:
                if seq.id == name_seq:
                    aa_at_pos = seq[pos - 1]
        except IndexError:
            print('The position {pos} is outside the sequence'.format(pos=pos))
            exit(1)
        return aa_at_pos

    def get_cat_in_range(self, pos1=None, pos2=None):
        """
        :param pos1: start position #from 1 until end
        :param pos2: end position #from 1 until end
        :return: list of the categories of amino acid at every position between pos1 and pos2
        """
        if pos1 is None:
            pos1 = 1
        if pos2 is None:
            pos2 = len(self.seqrefs[0])
        return self.list_of_categories[pos1 - 1:pos2]

    def get_aa_in_range(self, name_seq, pos1=None, pos2=None):
        """
        :param name_seq: name of the sequence
        :param pos1: beginning of the interval
        :param pos2: end of the interval
        :return: list of all the amino acids from the sequence in the interval of positions
        """
        if pos1 is None:
            pos1 = 1
        if pos2 is None:
            pos2 = len(self.seqeval[0])
        aa_in_range = []
        for pos in range(pos1, pos2 + 1):
            aa_in_range.append(set(self.get_aa_at_pos(pos, name_seq)))
        return aa_in_range

    def get_category_list(self, positions_list):
        """
        :param positions_list: list of the intervals of positions
        :return: list of the categories associated to the intervals of positions
        """
        list_of_categories = []
        for pos in range(len(positions_list)):
            if len(positions_list[pos]) == 1:
                list_of_categories.append([self.list_of_categories[pos - 1]])
            else:
                list_of_categories.append(
                    self.get_cat_in_range(positions_list[pos][0], positions_list[pos][1]))
        return list_of_categories

    def get_aa_list(self, positions_list, name_seq):
        """
        :param positions_list: list of the intervals of positions
        :param name_seq: name of the sequence
        :return: list of amino acids associated to the intervals of positions
        """
        list_of_aa = []
        for pos in positions_list:
            if len(pos) == 1:
                list_of_aa.append(
                    [set(self.get_aa_at_pos(pos[0], name_seq))])
            else:
                list_of_aa.append(
                    self.get_aa_in_range(name_seq, pos[0], pos[1]))
        return list_of_aa

    def entropy_pos_obs(self, pos):
        """
        :param pos: position in the alignment
        :return: the entropy associated to the position
        """
        return entropy(pk=[v for v in self.aa_ref_counts[pos - 1].values()], base=2)

    def entropy_background(self, method, gaps):
        """
        :param method: calculation method
        :param gaps: True if you want to consider gaps, False if not
        :return: the background entropy in the reference alignment
        """
        if method == 'database':
            ref_frq = {"A": 9.26, "Q": 3.75, "L": 9.91, "S": 6.63, "R": 5.80, "E": 6.16, "K": 4.88,
                       "T": 5.55, "N": 3.80, "G": 7.36, "M": 2.36, "W": 1.31, "D": 5.49, "H": 2.19,
                       "F": 3.91, "Y": 2.90, "C": 1.18, "I": 5.64, "P": 4.88, "V": 6.93}
            if gaps:
                ref_frq["-"] = 0
                for pos in self.aa_ref_counts:
                    ref_frq["-"] += pos["-"]
            pk = [v for v in ref_frq.values()]
        elif method == 'ref':
            count_all = {"A": 0, "R": 0, "N": 0, "D": 0, "C": 0, "E": 0, "Q": 0,
                         "G": 0, "H": 0, "I": 0, "L": 0, "K": 0, "M": 0, "F": 0, "P": 0, "S": 0,
                         "T": 0, "W": 0, "Y": 0, "V": 0}
            if gaps:
                count_all["-"] = 0
            for pos in self.aa_ref_counts:
                for i in count_all.keys():
                    count_all[i] += pos[i]
            pk = [v for v in count_all.values()]
        else:
            j = 21 if gaps else 20
            pk = [1 / j] * j
        return entropy(pk, base=2)

    def information_pos(self, pos, method, gaps):
        """
        :param pos: position in the alignment
        :param method: calculation method
        :param gaps: True if you want to consider gaps, False if not
        :return: the information at the position
        """
        return self.entropy_background(method, gaps) - self.entropy_pos_obs(pos)

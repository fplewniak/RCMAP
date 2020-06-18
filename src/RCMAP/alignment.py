from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from RCMAP.classification_aa import AAcategories
from scipy.stats import entropy
from RCMAP.utilities import get_weight, get_entropy_back


class Alignments:
    """
    Opens the file and processes the alignment of sequences.
    """

    def __init__(self, file, seqs_to_evaluate):
        self.file = file
        self.seqs_to_evaluate = seqs_to_evaluate
        try:
            self.alignment = AlignIO.read(file, "fasta")
        except FileNotFoundError:
            print('File {file} not found'.format(file=file))
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
        self.list_of_categories, self.list_of_cat_names, self.aa_observed = \
            self.determine_ref_categories()

    def count_aa_ref(self):
        """
        Counts the number of every amino acids at every position in the reference sequences, and
        lists them in a dictionary.

        :return: the count of amino acids at every position in all reference sequences
        """
        self.aa_ref_counts = [
            {"A": 0, "R": 0, "N": 0, "D": 0, "B": 0, "C": 0, "E": 0, "Q": 0, "Z": 0, "G": 0, "H": 0,
             "I": 0, "L": 0, "K": 0, "M": 0, "F": 0, "P": 0, "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0,
             "-": 0} for sub in range(self.alignment.get_alignment_length())]
        for pos in range(self.alignment.get_alignment_length()):
            for aa in self.aa_ref_counts[pos]:
                self.aa_ref_counts[pos][aa] = self.seqrefs[:, pos].count(aa)
        return self.aa_ref_counts

    def determine_ref_categories(self):
        """
        Recovers information about every position in the reference sequences.

        :return: the list of categories of amino acids at every position in the reference sequences,
         the list of the names of the categories and the dictionary of all the amino acids observed
         at every position and their count
        """
        self.list_of_categories = [set() for sub in range(self.alignment.get_alignment_length())]
        self.list_of_cat_names = [set() for sub in range(self.alignment.get_alignment_length())]
        self.aa_observed = [dict() for sub in range(self.alignment.get_alignment_length())]
        for pos in range(len(self.count_aa_ref())):
            self.list_of_categories[pos], self.list_of_cat_names[pos] = \
                AAcategories().find_category(
                    {aa for aa in self.aa_ref_counts[pos] if self.aa_ref_counts[pos][aa] > 0})
            for aa in self.aa_ref_counts[pos]:
                if self.aa_ref_counts[pos][aa] > 0:
                    self.aa_observed[pos][aa] = self.aa_ref_counts[pos][aa]
        return self.list_of_categories, self.list_of_cat_names, self.aa_observed


    def get_cat_name_at_pos(self, pos):
        """
        Recovers the name of the category of amino acids observed at a position in the reference
        sequences. Uses the position -1 because the user counts from 1 and python from 0.

        :param pos: position of the amino acid in seqrefs
        :return: the name of the category of amino acids observed in seqrefs
        """
        return self.list_of_cat_names[pos - 1]

    def get_cat_at_pos(self, pos, strict):
        """
        Recovers the category (set of amino acids) of amino acids observed at a position in the
        reference sequences. Uses the position -1 because the user counts from 1 and python
        from 0.

        :param strict: if True, the category is only the amino acids observed
        :param pos: position of the amino acid in seqrefs
        :return: the category (set) of amino acids observed in seqrefs
        """
        if strict:
            return {aa for aa in self.aa_observed[pos-1]}
        return self.list_of_categories[pos - 1]

    def get_aa_observed_at_pos(self, pos):
        """
        Recovers the amino acids observed at a position in the reference sequences and the count of
        these amino acids. The data are in a dictionary sorted from the amino acid the most present
        to the least present at this position.
        Uses the position -1 because the user counts from 1 and python from 0.

        :param pos:  position of the amino acids in seqrefs
        :return: all the amino acids observed in seqrefs at this position and the count of this
        amino acids in a sorted dictionary
        """
        return dict(sorted(self.aa_observed[pos-1].items(), key=lambda item: item[1], reverse=True))

    def get_aa_at_pos(self, pos, name_seq):
        """
        Recovers the amino acid present at this position in the given sequence.

        :param pos: position of the amino acid in seqref or seqeval
        :param name_seq: name of the sequence in seqref or seqeval
        :return: set() containing the amino acid at this position in the sequence
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

    def entropy_pos_obs(self, pos):
        """
        Calculates the entropy from the frequencies of the amino acids observed at a position.
        NB : the entropy function can take in parameter the accounts of amino acids or frequencies
        directly.

        :param pos: position in the alignment
        :return: the entropy associated to the position
        """
        return entropy(pk=[v for v in self.aa_ref_counts[pos - 1].values()], base=2)

    def entropy_background(self, method, gaps):
        """
        Calculates the background entropy from the frequencies of the amino acids given by the
        method. Can take into account the gaps in the frequencies if gaps = True.
        There is three methods :
        * database : the frequencies of the amino acids come from the bank UniprotKB,TrEMBL
          april 2020;
        * ref : the frequencies of the amino acids come from the average of the counts in the
          reference sequences;
        * equiprobable : the frequencies of the amino acids are all the same.
        NB : the entropy function can take in parameter the accounts of amino acids or frequencies
        directly.

        :param method: calculation method
        :param gaps: True if you want to consider gaps, False if not
        :return: the background entropy in the reference alignment
        """
        count_gaps = None
        if gaps:
            count_gaps = {'-': 0}
            for pos in self.aa_ref_counts:
                count_gaps["-"] += pos["-"]
        if method == 'ref':
            count_all = {"A": 0, "R": 0, "N": 0, "D": 0, "C": 0, "E": 0, "Q": 0,
                         "G": 0, "H": 0, "I": 0, "L": 0, "K": 0, "M": 0, "F": 0, "P": 0, "S": 0,
                         "T": 0, "W": 0, "Y": 0, "V": 0}
            for pos in self.aa_ref_counts:
                for i in count_all.keys():
                    count_all[i] += pos[i]
            if gaps:
                count_all.update(count_gaps)
            entropy_back = entropy([v for v in count_all.values()], base=2)
        else:
            entropy_back = get_entropy_back(method, count_gaps)
        return entropy_back

    def information_pos(self, pos, method, gaps, window, window_method):
        """
        Calculates the information carried by a position. The running window can be used to
        consider the environment of a position and smooth its information.

        :param window_method: calculation method of the weights at every position in the window
        :param window: number of positions to calculate the average of information, should be odd
        :param pos: position in the alignment
        :param method: calculation method of the frequencies in the background entropy
        :param gaps: True if you want to consider gaps, False if not
        :return: the information of the position
        """
        if window == 1:
            return self.entropy_background(method, gaps) - self.entropy_pos_obs(pos)
        weights = get_weight(window, window_method)
        w, info = int(window / 2), 0
        if pos - w < 1:
            k = - (pos - w) + 1  # To take into account the good weights associated with the
            # positions if max(pos - w, 1) = 1
        else:
            k = 0
        for i in range(max(pos - w, 1), min(pos + w, self.alignment.get_alignment_length()) + 1):
            info += (self.entropy_background(method, gaps) - self.entropy_pos_obs(i)) * weights[k]
            k += 1
        return info / sum(weights)

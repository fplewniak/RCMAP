import unittest
from Bio import AlignIO
from RCMAP.Classification_AA import AAcategories
from tests.read_alignment import Read_sequences


class MyTestCase(unittest.TestCase):

    def __init__(self,n=6,data="Test_reference_sequence.faa" ):
        """Save the number of reference sequences"""
        self.n= n
        self.data = data

    def test_reference_sequences(self):
        sq = AlignIO.read(self.data, "fasta")
        assert Read_sequences.reference_sequences() == sq

    def test_create_set_AA(self):
        assert Read_sequences.create_set_AA()== [{'M'}, {'G', 'S', 'D', 'H', 'P'}, {'D', '-', 'K'}, {'-', 'T'}, {'-', 'P'}, {'-', 'S'}, {'-', 'T'}, {'-', 'T'}, {'-', 'A'}, {'-', 'E'}, {'-', 'R'}, {'-', 'S'}, {'F', '-', 'H'}, {'-', 'L'}, {'F', 'I', '-', 'L'}, {'Q', '-', 'G', 'E'}]

    def test_create_set_categories(self):
        self.List_AA = [{'M'},{"-","K","D"},{"-","D"}]
        assert Read_sequences.create_set_categories() == [{'M'},{'E', 'K', 'H', 'R', 'D'},{"D"}]
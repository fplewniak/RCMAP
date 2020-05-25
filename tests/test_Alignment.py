import unittest
from Bio import AlignIO
from RCMAP.Classification_AA import AAcategories
from tests.read_alignment import Read_sequences


class MyTestCase(unittest.TestCase):

    def test_read_alignments(self):
        (seqref,seqeval) = RCMAP.Alignment().read_alignments("data/ArsM_aln.faa")
        assert [s.id for s in seqref] == ["gb|ACN39191.1","ref|XP_005539091.1","gb|PXF49941.1","gb|OSX76803.1","emb|SMG66974.1","ref|XP_005706547.1"]
        assert [s.id for s in seqeval] == ["Q968Z2","WP_045226361"]


#    def test_create_set_AA(self):
#        assert Read_sequences.create_set_AA()== [{'M'}, {'G', 'S', 'D', 'H', 'P'}, {'D', '-', 'K'}, {'-', 'T'}, {'-', 'P'}, {'-', 'S'}, {'-', 'T'}, {'-', 'T'}, {'-', 'A'}, {'-', 'E'}, {'-', 'R'}, {'-', 'S'}, {'F', '-', 'H'}, {'-', 'L'}, {'F', 'I', '-', 'L'}, {'Q', '-', 'G', 'E'}]

#    def test_create_set_categories(self):
#        self.List_AA = [{'M'},{"-","K","D"},{"-","D"}]
#        assert Read_sequences.create_set_categories() == [{'M'},{'E', 'K', 'H', 'R', 'D'},{"D"}]
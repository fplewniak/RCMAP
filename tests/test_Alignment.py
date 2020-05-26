import unittest
from RCMAP.Alignment import Alignments
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

class MyTestCase(unittest.TestCase):

    def test_read_alignments(self):
        (seqref, seqeval) = Alignments("ArsM_aln.faa",["WP_045226361.1", "Q969Z2"]).read_alignments()
        assert [s.id for s in seqref] == ["gb|ACN39191.1", "ref|XP_005539091.1", "gb|PXF49941.1", "gb|OSX76803.1",
                                          "emb|SMG66974.1", "ref|XP_005706547.1"]
        assert [s.id for s in seqeval] == ["WP_045226361.1", "Q969Z2"]

    def test_get_cat_at_pos(self):
        seqrefs = Alignments("ArsM_aln.faa",["WP_045226361.1", "Q969Z2"]).read_alignments()[0]
        assert Alignments("ArsM_aln.faa",["WP_045226361.1", "Q969Z2"]).get_cat_at_pos(3) == {'T'}

#    def test_create_set_categories(self):
#        self.List_AA = [{'M'},{"-","K","D"},{"-","D"}]
#        assert Read_sequences.create_set_categories() == [{'M'},{'E', 'K', 'H', 'R', 'D'},{"D"}]
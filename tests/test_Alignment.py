import unittest
from RCMAP.Alignment import Alignments
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

class MyTestCase(unittest.TestCase):

    def test_get_alignments(self):
        (seqref, seqeval) = Alignments("ArsM_aln.faa",["WP_045226361.1", "Q969Z2"]).get_alignments()
        assert [s.id for s in seqref] == ["gb|ACN39191.1", "ref|XP_005539091.1", "gb|PXF49941.1", "gb|OSX76803.1",
                                          "emb|SMG66974.1", "ref|XP_005706547.1"]
        assert [s.id for s in seqeval] == ["WP_045226361.1", "Q969Z2"]

    def test_get_cat_at_pos(self):
        assert Alignments("ArsM_aln.faa",["WP_045226361.1", "Q969Z2"]).get_cat_at_pos(4) == {'T'}

    def test_get_cat_in_range1(self):
        assert Alignments("ArsM_aln_part.faa",["WP_045226361.1", "Q969Z2"]).get_cat_in_range(1,3) == [{"M"},{"I","V","L","F","Y","W","H","M","K","T","G","A","C","P","S","N","D","E","Q","R"},{"D","E","K","R","H"}]
        assert Alignments("ArsM_aln_part.faa",["WP_045226361.1", "Q969Z2"]).get_cat_in_range(4,4)  == [{'T'}]

    def test_get_cat_in_range2(self):
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_cat_in_range(None, 2) == [{"M"},{"I","V","L","F","Y","W","H","M","K","T","G","A","C","P","S","N","D","E","Q","R"}]

    def test_get_cat_in_range3(self):
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_cat_in_range(4, None) == [{"T"},{"P"},{"S"},{"T"},{"T"},{"A"}]

    #def test_get_cat_in_range4(self):
    #    assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_cat_in_range(0, 10) == "Error"

    def test_get_positions_list(self):
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_positions_list(['3:10','8:25']) == [[3,10],[8,25]]

    def test_get_positions_list(self):
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_positions_list([':25' ,'32', '45:']) == [[None,25], [32], [45,None]]

    def test_get_category_list(self):
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_category_list([[1,2],[None,3]]) == [[{"M"},{"I","V","L","F","Y","W","H","M","K","T","G","A","C","P","S","N","D","E","Q","R"}],[{"M"},{"I","V","L","F","Y","W","H","M","K","T","G","A","C","P","S","N","D","E","Q","R"},{"D","E","K","R","H"}]]
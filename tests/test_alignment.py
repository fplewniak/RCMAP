import unittest
from RCMAP.alignment import Alignments
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from math import *


class MyTestCase(unittest.TestCase):

    def test_get_alignments(self):
        (seqref, seqeval) = Alignments("ArsM_aln.faa",
                                       ["WP_045226361.1", "Q969Z2"]).get_alignments()
        assert [s.id for s in seqref] == ["gb|ACN39191.1", "ref|XP_005539091.1", "gb|PXF49941.1",
                                          "gb|OSX76803.1",
                                          "emb|SMG66974.1", "ref|XP_005706547.1"]
        assert [s.id for s in seqeval] == ["WP_045226361.1", "Q969Z2"]

    def test_count_aa_ref(self):
        assert Alignments("ArsM_aln_part1.faa", ["WP_045226361.1", "Q969Z2"]).count_aa_ref() == [
            {"A": 0, "R": 0, "N": 0, "D": 0, "B": 0, "C": 0, "E": 0, "Q": 0, "Z": 0, "G": 0, "H": 0,
             "I": 0, "L": 0, "K": 0, "M": 6, "F": 0, "P": 0, "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0,
             "-": 0},
            {"A": 0, "R": 0, "N": 0, "D": 1, "B": 0, "C": 0, "E": 0, "Q": 0, "Z": 0, "G": 1, "H": 1,
             "I": 0, "L": 0, "K": 0, "M": 0, "F": 0, "P": 2, "S": 1, "T": 0, "W": 0, "Y": 0, "V": 0,
             "-": 0},
            {"A": 0, "R": 0, "N": 0, "D": 1, "B": 0, "C": 0, "E": 0, "Q": 0, "Z": 0, "G": 0, "H": 0,
             "I": 0, "L": 0, "K": 1, "M": 0, "F": 0, "P": 0, "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0,
             "-": 4}]

    def test_determine_ref_categories(self):
        assert Alignments("ArsM_aln_part1.faa",
                          ["WP_045226361.1", "Q969Z2"]).determine_ref_categories() == ([{'M'}, set(
            "IVLFYWHMKTGACPSNDEQR"), set("DEKRH")], [{'M'}, 'Any', 'Charged'], [{'M'},
                                                                                {'P', 'H', 'D', 'S',
                                                                                 'G'},
                                                                                {'-', 'K', 'D'}])

    def test_get_aa_at_pos(self):
        assert Alignments("ArsM_aln.faa", ["WP_045226361.1", "Q969Z2"]).get_aa_at_pos(37,
                                                                                      "WP_045226361.1") == 'M'
        assert Alignments("ArsM_aln.faa", ["WP_045226361.1", "Q969Z2"]).get_aa_at_pos(16,
                                                                                      "emb|SMG66974.1") == 'G'

    #    def test_get_cat_at_pos(self):
    #        assert Alignments("ArsM_aln.faa",["WP_045226361.1", "Q969Z2"]).get_cat_at_pos(4) == {'T'}

    def test_get_cat_in_range(self):
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_cat_in_range(1,
                                                                                              3) == [
                   {"M"},
                   {"I", "V", "L", "F", "Y", "W", "H", "M", "K", "T", "G", "A", "C", "P", "S", "N",
                    "D", "E", "Q", "R"}, {"D", "E", "K", "R", "H"}]
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_cat_in_range(4,
                                                                                              4) == [
                   {'T'}]
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_cat_in_range(None,
                                                                                              2) == [
                   {"M"},
                   {"I", "V", "L", "F", "Y", "W", "H", "M", "K", "T", "G", "A", "C", "P", "S", "N",
                    "D", "E", "Q", "R"}]
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_cat_in_range(4,
                                                                                              None) == [
                   {"T"}, {"P"}, {"S"}, {"T"}, {"T"}, {"A"}]

    # def test_get_cat_in_range4(self):
    #    assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_cat_in_range(0, 10) == "Error"

    def test_get_aa_in_range(self):
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_aa_in_range(
            "Q969Z2", 1, 3) == [{"M"}, {"A"}, {"-"}]
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_aa_in_range(
            "Q969Z2", None, 2) == [{"M"}, {"A"}]
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_aa_in_range(
            "Q969Z2", 4,
            None) == [{"-"}, {"-"}, {"-"}, {"-"}, {"-"}, {"-"}]

    def test_get_category_list(self):
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_category_list(
            [[1, 2], [None, 3]]) == [[{"M"},
                                      {"I", "V", "L", "F", "Y", "W", "H", "M", "K", "T", "G", "A",
                                       "C", "P", "S", "N", "D", "E", "Q", "R"}], [{"M"},
                                                                                  {"I", "V", "L",
                                                                                   "F", "Y", "W",
                                                                                   "H", "M", "K",
                                                                                   "T", "G", "A",
                                                                                   "C", "P", "S",
                                                                                   "N", "D", "E",
                                                                                   "Q", "R"},
                                                                                  {"D", "E", "K",
                                                                                   "R", "H"}]]
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_category_list(
            [[1]]) == [[{"A"}]]

    def test_get_aa_list(self):
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_aa_list(
            [[1, 2], [None, 3]], "Q969Z2") == [[{'M'}, {'A'}], [{'M'}, {'A'}, {'-'}]]

    def test_entropy_pos_obs(self):
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).entropy_pos_obs(1) == 0
        assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).entropy_pos_obs(
            3) == - ((1 / 6 * log(1 / 6, 2)) * 2 + 4 / 6 * log(4 / 6, 2))

    def test_entropy_pos_background(self):
        assert Alignments("ArsM_aln_part.faa",
                          ["WP_045226361.1", "Q969Z2"]).entropy_pos_background() == - (
                9.26 * log(9.26, 2) + 3.75 * log(3.75, 2) + 9.91 * log(9.91, 2) + 6.63 * log(6.63,
                                                                                             2) + 5.8 * log(
            5.8, 2) + 6.16 * log(6.16, 2) + 4.88 * log(4.88, 2) + 5.55 * log(5.55, 2) + 3.8 * log(
            3.8, 2) + 7.36 * log(7.36, 2) + 2.36 * log(2.36, 2) + 1.31 * log(1.31, 2) + 5.49 * log(
            5.49, 2) + 2.19 * log(2.19, 2) + 3.91 * log(3.91, 2) + 2.9 * log(2.9, 2) + 1.18 * log(
            1.18, 2) + 5.64 * log(5.64, 2) + 4.88 * log(4.88, 2) + 6.93 * log(6.93, 2))
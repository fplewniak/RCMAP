from RCMAP.alignment import Alignments
from math import *
from scipy.stats import entropy


def test_get_alignments():
    """
    test the recovery of the reference sequences and evaluates sequences
    """
    (seqref, seqeval) = Alignments("ArsM_aln.faa",
                                   ["WP_045226361.1", "Q969Z2"]).get_alignments()
    assert [s.id for s in seqref] == ["gb|ACN39191.1", "ref|XP_005539091.1", "gb|PXF49941.1",
                                      "gb|OSX76803.1",
                                      "emb|SMG66974.1", "ref|XP_005706547.1"]
    assert [s.id for s in seqeval] == ["WP_045226361.1", "Q969Z2"]


def test_count_aa_ref():
    """
    test the count of the amino acids at all positions
    """
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


def test_determine_ref_categories():
    """
    test the list of categories, the list of categories sets, the set of amino acids and
     the count of every amino acid in all positions of the reference sequences
    """
    assert Alignments("ArsM_aln_part1.faa",
                      ["WP_045226361.1", "Q969Z2"]).determine_ref_categories() == ([{'M'}, set(
        "IVLFYWHMKTGACPSNDEQR"), set("DEKRH")], [{'M'}, 'Any', 'Charged'], [['M'],
                                                                            ['D', 'G', 'H', 'P',
                                                                             'S'],
                                                                            ['D', 'K', '-']],
                                                                                   [[6],
                                                                                    [1, 1, 1, 2,
                                                                                     1],
                                                                                    [1, 1, 4]],
                                                                                   [{'M'},
                                                                                    {'S', 'P', 'G',
                                                                                     'H', 'D'},
                                                                                    {'-', 'K',
                                                                                     'D'}])


def test_get_aa_observed_at_pos():
    """
    """
    assert Alignments("ArsM_aln.faa", ["WP_045226361.1", "Q969Z2"]).get_aa_observed_at_pos(1) == {
        'M': 6}
    assert Alignments("ArsM_aln.faa", ["WP_045226361.1", "Q969Z2"]).get_aa_observed_at_pos(3) == {
        '-': 4, 'D': 1, 'K': 1}


def test_get_aa_at_pos():
    """
    test the amino acid at a position in the sequence
    """
    assert Alignments("ArsM_aln.faa", ["WP_045226361.1", "Q969Z2"]).get_aa_at_pos(37,
                                                                                  "WP_045226361.1") \
           == 'M'
    assert Alignments("ArsM_aln.faa", ["WP_045226361.1", "Q969Z2"]).get_aa_at_pos(16,
                                                                                  "emb|SMG66974.1") \
           == 'G'



def test_get_cat_in_range():
    """
    test the categories of the amino acids in the reference sequences in the given interval
    """
    assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_cat_in_range(1,
                                                                                          3) == [
               {"M"},
               {"I", "V", "L", "F", "Y", "W", "H", "M", "K", "T", "G", "A", "C", "P", "S",
                "N",
                "D", "E", "Q", "R"}, {"D", "E", "K", "R", "H"}]
    assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_cat_in_range(4,
                                                                                          4) == [
               {'T'}]


def test_get_aa_in_range():
    """
    test the amino acids in the given sequence in the given interval
    """
    assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_aa_in_range(
        "Q969Z2", 1, 3) == [{"M"}, {"A"}, {"-"}]


def test_get_category_list():
    """
    test the list of the categories in the reference sequences in the intervals given
    """
    assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_category_list(
        [[1, 2], [1, 3]]) == [[{"M"},
                                  {"I", "V", "L", "F", "Y", "W", "H", "M", "K", "T", "G",
                                   "A",
                                   "C", "P", "S", "N", "D", "E", "Q", "R"}], [{"M"},
                                                                              {"I", "V",
                                                                               "L",
                                                                               "F", "Y",
                                                                               "W",
                                                                               "H", "M",
                                                                               "K",
                                                                               "T", "G",
                                                                               "A",
                                                                               "C", "P",
                                                                               "S",
                                                                               "N", "D",
                                                                               "E",
                                                                               "Q", "R"},
                                                                              {"D", "E",
                                                                               "K",
                                                                               "R", "H"}]]
    assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_category_list(
        [[1]]) == [[{"A"}]]


def test_get_aa_list():
    """
    test list of the amino acids in the given sequence in the interval given
    """
    assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).get_aa_list(
        [[1, 2], [1, 3]], "Q969Z2") == [[{'M'}, {'A'}], [{'M'}, {'A'}, {'-'}]]


def test_entropy_pos_obs():
    """
    test the entropy calculated at a position
    """
    assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).entropy_pos_obs(
        1) == 0
    assert Alignments("ArsM_aln_part.faa", ["WP_045226361.1", "Q969Z2"]).entropy_pos_obs(
        3) == - ((1 / 6 * log(1 / 6, 2)) * 2 + 4 / 6 * log(4 / 6, 2))


def test_entropy_background():
    """
    test the calculation of the background entropy with different methods
    """
    assert Alignments("ArsM_aln_part.faa",
                      ["WP_045226361.1", "Q969Z2"]).entropy_background('database',
                                                                       True) == entropy(
        pk=[9.26, 3.75, 9.91, 6.63, 5.80, 6.16, 4.88,
            5.55, 3.80, 7.36, 2.36, 1.31, 5.49, 2.19,
            3.91, 2.90, 1.18, 5.64, 4.88, 6.93, 34], qk=None, base=2)
    assert Alignments("ArsM_aln_part.faa",
                      ["WP_045226361.1", "Q969Z2"]).entropy_background('database',
                                                                       False) == entropy(
        pk=[9.26, 3.75, 9.91, 6.63, 5.80, 6.16, 4.88,
            5.55, 3.80, 7.36, 2.36, 1.31, 5.49, 2.19,
            3.91, 2.90, 1.18, 5.64, 4.88, 6.93], qk=None, base=2)
    assert round(Alignments("ArsM_aln_part.faa",
                            ["WP_045226361.1", "Q969Z2"]).entropy_background('equiprobable',
                                                                             False),
                 6) == round(- (
            1 / 20 * log2(1 / 20)) * 20, 6)
    assert round(Alignments("ArsM_aln_part1.faa",
                            ["WP_045226361.1", "Q969Z2"]).entropy_background('equiprobable',
                                                                             True),
                 6) == round(- (
            1 / 21 * log2(1 / 21)) * 21, 6)
    assert round(Alignments("ArsM_aln_part1.faa",
                            ["WP_045226361.1", "Q969Z2"]).entropy_background('ref', False),
                 6) == round(- (
            (2 / 14 * log2(2 / 14)) * 2 + (1 / 14 * log2(1 / 14)) * 4 + 6 / 14 * log2(
        6 / 14)), 6)
    assert round(Alignments("ArsM_aln_part1.faa",
                            ["WP_045226361.1", "Q969Z2"]).entropy_background('ref',
                                                                             True),
                 6) == round(- (
            (2 / 18 * log2(2 / 18)) * 2 + (1 / 18 * log2(1 / 18)) * 4 + 4 / 18 * log2(
        4 / 18) + 6 / 18 * log2(6 / 18)), 6)

def test_information_pos():
    assert round(Alignments("ArsM_aln_part1.faa",
                            ["WP_045226361.1", "Q969Z2"]).information_pos(2, 'database', False, 5), 2) == 2.99
    assert round(Alignments("ArsM_aln_part1.faa",
                            ["WP_045226361.1", "Q969Z2"]).information_pos(2, 'database', False, 1), 2) == 1.9

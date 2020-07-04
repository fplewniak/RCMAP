import pytest
from math import *
from RCMAP.utilities import get_positions_list, compatibility, sort_categories, summary_info, \
    get_weight, get_entropy_back
from scipy.stats import entropy


def test_get_positions_list():
    """
    Tests that positions like : ['3:10', '8:25' ,'32', '45:'] are transformed in positions like :
    [[3, 10], [8, 25], [32, 32], [45, pos_max]].
    """
    assert get_positions_list(['3:10', '8:25'], 100) == [[3, 10], [8, 25]]
    assert get_positions_list([':25', '32', '45:'], 100) == [[1, 25], [32, 32], [45, 100]]
    assert get_positions_list(['20:10'], 100) == [[10, 20]]
    with pytest.raises(SystemExit):
        assert get_positions_list(['1:-5'], 100)
    with pytest.raises(SystemExit):
        assert get_positions_list(['-10'], 100)


def test_compatibility():
    """
    Tests the compatibility of different amino acid with categories.
    """
    assert compatibility({"B"}, set("ASCGPNDTV"), False) is True
    assert compatibility({"Z"}, set("DEKRHQNSCTYW"), False) is True
    assert compatibility({"-"}, set("IVL"), False) is False
    assert compatibility({"-"}, set("IVL"), True) is True
    assert compatibility({"-"}, set("-VL"), True) is True
    assert compatibility({"G"}, set("RKH"), True) is False
    assert compatibility({"G"}, set("ASCGPNDTV"), False) is True
    assert compatibility({"G"}, set("RKH"), False) is False


def test_sort_categories():
    """
    Tests that the sorted list of categories is the good one.
    """
    assert sort_categories() == ['Negative', 'Positive', 'Aliphatic', 'Tiny',
                                 'Aromatic', 'Charged', 'Small', 'Polar',
                                 'Hydrophobic']


def test_summary_info():
    """
    Tests the number of true in the list and the information associated (same for false).
    """
    assert summary_info([False, True, True, False], [2, 1, 5, 3]) == (2, 6, 2, 5)


def test_get_weight():
    """
    Tests the list of weights associated to the different methods.
    """
    assert [w for w in get_weight(5, 'flat')] == [1., 1., 1., 1., 1.]
    assert [w for w in get_weight(3, 'bartlett')] == [
        2 / (3 - 1) * ((3 - 1) / 2 - abs(0 - (3 - 1) / 2)),
        2 / (3 - 1) * ((3 - 1) / 2 - abs(1 - (3 - 1) / 2)),
        2 / (3 - 1) * ((3 - 1) / 2 - abs(0 - (3 - 1) / 2))]
    assert [w for w in get_weight(3, 'hamming')] == [0.54 - 0.46 * cos(2 * pi * 0 / (3 - 1)),
                                                     0.54 - 0.46 * cos(2 * pi * 1 / (3 - 1)),
                                                     0.54 - 0.46 * cos(2 * pi * 2 / (3 - 1))]
    assert [w for w in get_weight(3, 'hanning')] == [0.5 - 0.5 * cos(2 * pi * 0 / (3 - 1)),
                                                     0.5 - 0.5 * cos(2 * pi * 1 / (3 - 1)),
                                                     0.5 - 0.5 * cos(2 * pi * 2 / (3 - 1))]


def test_get_entropy_back():
    """
    Tests the calculation of the background entropy with different methods.
    """
    assert get_entropy_back('equiprobable', None) == entropy([1 / 20] * 20, base=2)
    assert get_entropy_back('equiprobable', 14) == entropy([1 / 21] * 21, base=2)
    assert get_entropy_back('database', None) == entropy(
        [v for v in {"A": 9.26, "Q": 3.75, "L": 9.91, "S": 6.63, "R": 5.80, "E": 6.16, "K": 4.88,
                     "T": 5.55, "N": 3.80, "G": 7.36, "M": 2.36, "W": 1.31, "D": 5.49, "H": 2.19,
                     "F": 3.91, "Y": 2.90, "C": 1.18, "I": 5.64, "P": 4.88, "V": 6.93}.values()],
        base=2)
    assert get_entropy_back('database', {'-': 20}) == entropy(
        [v for v in {"A": 9.26, "Q": 3.75, "L": 9.91, "S": 6.63, "R": 5.80, "E": 6.16, "K": 4.88,
                     "T": 5.55, "N": 3.80, "G": 7.36, "M": 2.36, "W": 1.31, "D": 5.49, "H": 2.19,
                     "F": 3.91, "Y": 2.90, "C": 1.18, "I": 5.64, "P": 4.88, "V": 6.93,
                     "-": 20}.values()],
        base=2)

import pytest
from math import *
from RCMAP.utilities import get_positions_list, compatibility, sort_categories, summary_info, \
    get_weight


def test_get_positions_list():
    """
    test : positions like ['3:10', '8:25' ,'32', '45:'] -> positions like [[3, 10], [8, 25],
    [32, 32], [45, pos_max]]
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
    tests the compatibility of different amino acid with categories
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
    tests that the sorted list of categories is the good one
    """
    assert sort_categories() == ['Negative', 'Positive', 'Aliphatic', 'Tiny',
                                 'Aromatic', 'Charged', 'Small', 'Polar',
                                 'Hydrophobic']


def test_summary_info():
    """
    tests the number of true in the list and the information associated (same for false)
    """
    assert summary_info([False, True, True, False], [2, 1, 5, 3]) == (2, 6, 2, 5)


def test_get_weight():
    """
    tests the list of weights associated to the different methods
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

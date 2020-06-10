from RCMAP.utilities import get_positions_list, compatibility, sort_categories, summary_info


def test_get_positions_list():
    """
    test : positions like ['3:10', '8:25' ,'32', '45:'] -> positions like [[3,10], [8,25], [32],
    [45,None]]
    """
    assert get_positions_list(['3:10', '8:25']) == [[3, 10], [8, 25]]
    assert get_positions_list([':25', '32', '45:']) == [[None, 25], [32, 32], [45, None]]


def test_sort_categories():
    """
    test that the sorted list of categories is the good one
    """
    assert sort_categories() == ['Negative', 'Positive', 'Aliphatic', 'Tiny',
                                 'Aromatic', 'Charged', 'Small', 'Polar',
                                 'Hydrophobic']


def test_summary_info():
    """
    test the number of true in the list and the information associated (same for false)
    """
    assert summary_info([False, True, True, False], [2, 1, 5, 3]) == (2, 6, 2, 5)

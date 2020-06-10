from RCMAP.classification_aa import AAcategories


def test_find_category():
    """
    test the category of a set of amino acids
    """
    assert AAcategories().find_category({"A"}) == ({"A"}, {"A"})
    assert AAcategories().find_category({}) == ({}, {})
    assert AAcategories().find_category({"I", "V", "A", "R"}) == (
            set("IVLFYWHMKTGACPSNDEQR"), 'Any')
    assert AAcategories().find_category({"-", "-"}) == (set(), set())
    assert AAcategories().find_category({"R", "H"}) == ({'R', 'H', 'K'}, 'Positive')
    assert AAcategories().find_category({"B"}) == (
            {'P', 'S', 'N', 'D', 'T', 'V', 'G', 'C', 'A'}, 'Small')
    assert AAcategories().find_category({"Z"}) == (
            {'Q', 'S', 'D', 'T', 'W', 'R', 'H', 'K', 'Y', 'C', 'N', 'E'}, 'Polar')



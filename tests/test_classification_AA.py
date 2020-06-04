import unittest
import pytest

from RCMAP.classification_aa import AAcategories

class MyTestCase(unittest.TestCase):

    def test_sort_categories(self):
        assert AAcategories().sort_categories() == ['Negative', 'Positive', 'Aliphatic', 'Tiny', 'Aromatic', 'Charged', 'Small', 'Polar', 'Hydrophobic']


    def test_find_category1(self):
        assert AAcategories().find_category({"A"}) == ({"A"}, {"A"})

    def test_find_category2(self):
        assert AAcategories().find_category({}) == ({}, {})

    def test_find_category3(self):
        assert AAcategories().find_category({"I","V","A","R"}) == (set("IVLFYWHMKTGACPSNDEQR"), 'Any')

    def test_find_category4(self):
        assert AAcategories().find_category({"-","-"}) == (set(), set())

    def test_find_category5(self):
        assert AAcategories().find_category({"R","H"}) == ({'R', 'H', 'K'}, 'Positive')

    def test_find_category6(self):
        assert AAcategories().find_category({"B"}) == ({'P', 'S', 'N', 'D', 'T', 'V', 'G', 'C', 'A'}, 'Small')

    def test_find_category7(self):
        assert AAcategories().find_category({"Z"}) == ({'Q', 'S', 'D', 'T', 'W', 'R', 'H', 'K', 'Y', 'C', 'N', 'E'}, 'Polar')


    def test_compatibility(self):
        assert AAcategories().compatibility({"B"},set("ASCGPNDTV")) == True

    def test_compatibility1(self):
        assert AAcategories().compatibility({"Z"},set("DEKRHQNSCTYW")) == True

    def test_compatibility2(self):
        assert AAcategories().compatibility({"-"}, set("IVL")) == False

    def test_compatibility3(self):
         assert AAcategories().compatibility({"G"}, set("ASCGPNDTV")) == True

    def test_compatibility3(self):
         assert AAcategories().compatibility({"G"}, set("RKH")) == False
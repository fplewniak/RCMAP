import unittest
import pytest

from RCMAP.Classification_AA import AAcategories

class MyTestCase(unittest.TestCase):

    def test_sort_categories(self):
        assert AAcategories().sort_categories() == ['Negative', 'Positive', 'Aliphatic', 'Tiny', 'Aromatic', 'Charged', 'Small', 'Polar', 'Hydrophobic']


    def test_find_category1(self):
        assert AAcategories().find_category({"A"}) == {"A"}

    def test_find_category2(self):
        assert AAcategories().find_category({}) == {}

    def test_find_category3(self):
        assert AAcategories().find_category({"I","V","A","R"}) == set("IVLFYWHMKTGACPSNDEQR")

    def test_find_category4(self):
        assert AAcategories().find_category({"-","-"}) == set()

    def test_find_category5(self):
        assert AAcategories().find_category({"R","H"}) == {"R","K","H"}

    def test_find_category6(self):
        assert AAcategories().find_category({"B"}) == set("ASCGPNDTV")

    def test_find_category7(self):
        assert AAcategories().find_category({"Z"}) == set("DEKRHQNSCTYW")


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
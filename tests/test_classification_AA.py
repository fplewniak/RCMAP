import unittest
import pytest

from RCMAP.Classification_AA import AAcategories

class MyTestCase(unittest.TestCase):

    def test_sort_categories(self):
        assert AAcategories().sort_categories() == ['Negative', 'Positive', 'Aliphatic', 'Tiny', 'Aromatic', 'Charged', 'Small', 'Polar', 'Hydrophobic']


    def test_find_category1(self):
        ens = {"A"}
        assert AAcategories().find_category(ens) == {"A"}

    def test_find_category2(self):
        ens = {}
        assert AAcategories().find_category(ens) == {}

    def test_find_category3(self):
        ens = {"I","V","A","R"}
        assert AAcategories().find_category(ens) == set("IVLFYWHMKTGACPSNDEQR")

    def test_find_category4(self):
        ens = {"-","-"}
        assert AAcategories().find_category(ens) == set()

    def test_find_category5(self):
        ens = {"R","H"}
        assert AAcategories().find_category(ens) == {"R","K","H"}

    def test_find_category6(self):
        ens = {"B"}
        assert AAcategories().find_category(ens) == set("ASCGPNDTV")

    def test_find_category7(self):
        ens = {"Z"}
        assert AAcategories().find_category(ens) == set("DEKRHQNSCTYW")


    def test_compatibility(self):
        AA= {"B"}
        category = set("ASCGPNDTV")
        assert AAcategories().compatibility(AA,category) == True

    def test_compatibility1(self):
        AA= {"Z"}
        category = set("DEKRHQNSCTYW")
        assert AAcategories().compatibility(AA,category) == True

    def test_compatibility2(self):
        AA = {"-"}
        category = set("IVL")
        assert AAcategories().compatibility(AA, category) == False

    def test_compatibility3(self):
         AA = {"G"}
         category = set("ASCGPNDTV")
         assert AAcategories().compatibility(AA, category) == True

    def test_compatibility3(self):
         AA = {"G"}
         category = set("RKH")
         assert AAcategories().compatibility(AA, category) == False
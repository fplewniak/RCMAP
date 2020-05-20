import unittest
import pytest

from RCMAP.Classification_AA import AAcategories

class MyTestCase(unittest.TestCase):

    def test_sort_categories():
        assert AAcategories.sort_categories(self) == ['Negative', 'Positive', 'Aliphatic', 'Tiny', 'Aromatic', 'Charged', 'Small', 'Polar', 'Hydrophobic']


objet = MyTestCase()
print(objet.test_sort_categories())

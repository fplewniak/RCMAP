import unittest
import pytest

from RCMAP.Classification_AA import AAcategories

class MyTestCase(unittest.TestCase):

    def test_sort_categories(self):
        assert AAcategories().sort_categories() == ['Negative', 'Positive', 'Aliphatic', 'Tiny', 'Aromatic', 'Charged', 'Small', 'Polar', 'Hydrophobic']

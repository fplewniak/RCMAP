import unittest
from RCMAP.information_measurement import *

class MyTestCase(unittest.TestCase):

    def test_frequence_aa(self):
        assert frequence_aa('M',1) == 1.0



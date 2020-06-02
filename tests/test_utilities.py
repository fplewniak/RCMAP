import unittest
from RCMAP.Utilities import get_positions_list

class MyTestCase(unittest.TestCase):

    def test_get_positions_list(self):
        assert get_positions_list(['3:10','8:25']) == [[3,10],[8,25]]

    def test_get_positions_list(self):
        assert get_positions_list([':25' ,'32', '45:']) == [[None,25], [32,32], [45,None]]

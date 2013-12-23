#!/usr/bin/env python

import unittest
from src.translate import *

class TestTranslate(unittest.TestCase):

    def setUp(self):
        foo = 'foo'

    def test_valid_seq(self):
        self.assertTrue(valid_seq('actg'))
        self.assertFalse(valid_seq('acth'))
        self.assertFalse(valid_seq('ac'))

    def test_has_start_codon(self):
        self.assertTrue(has_start_codon('auggattaca'))
        self.assertFalse(has_start_codon('guggattaca')) # currently no support for alternate start codons

    def test_has_stop_codon(self):
        self.assertTrue(has_stop_codon('gattacatag'))
        self.assertTrue(has_stop_codon('gattacataa'))
        self.assertTrue(has_stop_codon('gattacatga'))
        self.assertFalse(has_stop_codon('gattacaact'))
        
        


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranslate))
    return suite

if __name__ == '__main__':
    unittest.main()

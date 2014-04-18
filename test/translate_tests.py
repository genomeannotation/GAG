#!/usr/bin/env python

import unittest
from src.translate import *

class TestTranslate(unittest.TestCase):

    def test_translate_one_codon(self):
        self.assertEquals('T', translate('act', '+', 1))
        self.assertEquals('T', translate('acta', '+', 1))
        self.assertEquals('T', translate('actac', '+', 1))
        self.assertEquals('T', translate('tact', '+', 2))
        self.assertEquals('T', translate('ctact', '+', 3))

        self.assertEquals('S', translate('act', '-', 1))
        self.assertEquals('S', translate('tact', '-', 1))
        self.assertEquals('S', translate('ctact', '-', 1))
        self.assertEquals('S', translate('acta', '-', 2))
        self.assertEquals('S', translate('actac', '-', 3))

    def test_translate_error_handling(self):
        sys.stderr.write("test should generate an error message:\n")
        self.assertFalse(translate('ac', 1, '+'))
        sys.stderr.write("test should generate an error message:\n")
        self.assertFalse(translate('acx', 1, '+'))
        sys.stderr.write("test should generate an error message:\n")
        self.assertFalse(translate('act', 0, '+'))
        sys.stderr.write("test should generate an error message:\n")
        self.assertFalse(translate('act', 1, 'u'))

    def test_translate_longer_sequence(self):
        test_seq = 'CATGACAGAAGATATTTC'
        self.assertEquals('HDRRYF', translate(test_seq, '+', 1,))
        self.assertEquals('MTEDI', translate(test_seq, '+', 2))
        self.assertEquals('*QKIF', translate(test_seq, '+', 3))
        self.assertEquals('EISSVM', translate(test_seq, '-', 1))
        self.assertEquals('KYLLS', translate(test_seq, '-', 2))
        self.assertEquals('NIFCH', translate(test_seq, '-', 3))
        
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

    def test_reverse_complement(self):
        self.assertEquals('C', reverse_complement('G'))
        self.assertEquals('CAT', reverse_complement('ATG'))
        
        
##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranslate))
    return suite

if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python

import unittest
from fasta import Fasta

class TestFasta(unittest.TestCase):


    def test_fasta(self):
        fasta = Fasta()
        fasta.readString('>seq1\nATGCCGTA\n>seq2\nAGGTCC\n>seq3\nGGGGGG')
        fasta.trimSeq('seq2', 2, 3)
        self.assertEqual(fasta.writeString(), '>seq1\nATGCCGTA\n>seq2\nATCC\n>seq3\nGGGGGG')
        fasta.trimSeq('seq2', 1, 4)
        self.assertEqual(fasta.writeString(), '>seq1\nATGCCGTA\n>seq3\nGGGGGG')


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFasta))
    return suite

if __name__ == '__main__':
    unittest.main()

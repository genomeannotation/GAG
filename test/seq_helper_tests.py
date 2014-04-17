#!/usr/bin/env python

import unittest
from mock import Mock
from src.seq_helper import SeqHelper

class TestSeqHelper(unittest.TestCase):

    def setUp(self):
        self.helper = SeqHelper("nnnGATTACAnnnnnNNNNNgattacaNNN")

    def test_initialize(self):
        self.assertEquals("nnnGATTACAnnnnnNNNNNgattacaNNN", self.helper.full_sequence)

    def test_mrna_to_fasta(self):
        mrna = Mock()
        mrna.identifier = "foo_mrna"
        mrna.exon = Mock()
        mrna.exon.indices = [[4, 10], [21, 27]]
        mrna.strand = '+'
        expected = "foo_mrna\nGATTACAgattaca\n"
        self.assertEquals(expected, self.helper.mrna_to_fasta(mrna))


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSeqHelper))
    return suite

if __name__ == '__main__':
    unittest.main()

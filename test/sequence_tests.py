#!/usr/bin/env python

import unittest
from mock import Mock
from src.sequence import Sequence

class TestSequence(unittest.TestCase):

    def setUp(self):
        self.seq1 = Sequence("seq1", "GATTACA")

    def test_string(self):
        expected = "Sequence seq1 of length 7 containing 0 genes\n"
        self.assertEquals(expected, str(self.seq1))

    def test_how_many_Ns_forward(self):
        badseq = Sequence('seq1', 'NNnNNGATTACA')
        self.assertEqual(5, badseq.how_many_Ns_forward(1))

    def test_how_many_Ns_forward_returns_zero_if_no_Ns(self):
        badseq = Sequence('seq2', 'GATTACA')
        self.assertEqual(0, badseq.how_many_Ns_forward(3))

    def test_how_many_Ns_backward(self):
        badseq = Sequence('seq3', 'gattaNnN')
        self.assertEqual(3, badseq.how_many_Ns_backward(8))

    def test_how_many_Ns_backward_returns_zero_if_no_Ns(self):
        self.assertEqual(0, self.seq1.how_many_Ns_backward(3))

    def test_add_gene(self):
        self.assertEqual(0, len(self.seq1.genes))
        mockgene = Mock()
        self.seq1.add_gene(mockgene)
        self.assertEqual(1, len(self.seq1.genes))

    def test_remove_gene(self):
        mockgene = Mock()
        mockgene.identifier = "foo_gene"
        self.seq1.add_gene(mockgene)
        self.assertEqual(1, len(self.seq1.genes))
        self.seq1.remove_gene("foo_gene")
        self.assertEqual(0, len(self.seq1.genes))

    def test_remove_gene_fails_if_no_match(self):
        mockgene = Mock()
        mockgene.identifier = "foo_gene"
        self.seq1.add_gene(mockgene)
        self.assertEqual(1, len(self.seq1.genes))
        self.seq1.remove_gene("bar_gene")
        self.assertEqual(1, len(self.seq1.genes))

    def test_trim(self):
        self.assertEquals("GATTACA", self.seq1.bases)
        self.seq1.trim(1, 4)
        self.assertEquals("ACA", self.seq1.bases)

    def test_trim_trims_gene(self):
        mockgene = Mock()
        mockgene.indices = [3, 7]
        self.seq1.add_gene(mockgene)
        self.assertEquals(1, len(self.seq1.genes))
        self.seq1.trim(1, 4)
        self.assertEquals(0, len(self.seq1.genes))

    def test_get_subseq(self):
        pass

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSequence))
    return suite

if __name__ == '__main__':
    unittest.main()

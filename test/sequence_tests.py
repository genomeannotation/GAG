#!/usr/bin/env python

import unittest
from mock import Mock
from src.sequence import Sequence

class TestSequence(unittest.TestCase):

    def setUp(self):
        self.seq1 = Sequence("seq1", "GATTACA")

    def add_mock_gene(self):
        mockgene = Mock()
        mockgene.identifier = "foo_gene"
        self.seq1.add_gene(mockgene)
        
    def add_mock_gene_with_1_mrna(self):
        mockgene = Mock()
        mockgene.identifier = "foo_gene"
        mockgene.mrnas = [Mock()]
        mockgene.mrnas[0].cds = Mock()
        mockgene.mrnas[0].exon = Mock()
        self.seq1.add_gene(mockgene)
        
    def add_mock_gene_with_2_mrnas(self):
        mockgene = Mock()
        mockgene.identifier = "foo_gene"
        mockgene.mrnas = [Mock(), Mock()]
        mockgene.mrnas[0].cds = None
        mockgene.mrnas[0].exon = None
        mockgene.mrnas[1].cds = Mock()
        mockgene.mrnas[1].exon = Mock()
        self.seq1.add_gene(mockgene)

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
        self.add_mock_gene()
        self.assertEqual(1, len(self.seq1.genes))

    def test_remove_gene(self):
        self.add_mock_gene()
        self.assertEqual(1, len(self.seq1.genes))
        self.seq1.remove_gene("foo_gene")
        self.assertEqual(0, len(self.seq1.genes))

    def test_remove_gene_fails_if_no_match(self):
        self.add_mock_gene()
        self.assertEqual(1, len(self.seq1.genes))
        self.seq1.remove_gene("bar_gene")
        self.assertEqual(1, len(self.seq1.genes))

    def test_trim_region(self):
        self.assertEquals("GATTACA", self.seq1.bases)
        self.seq1.trim_region(1, 4)
        self.assertEquals("ACA", self.seq1.bases)

    def test_trim_region_trims_gene(self):
        # TODO don't remove gene, trim its indices and verify?
        self.add_mock_gene()
        self.seq1.genes[0].indices = [3, 7]
        self.assertEquals(1, len(self.seq1.genes))
        self.seq1.trim_region(1, 4)
        self.assertEquals(0, len(self.seq1.genes))

    def test_get_subseq(self):
        self.assertEquals("ATTA", self.seq1.get_subseq(2, 5))

    def test_to_tbl(self):
        self.add_mock_gene()
        self.seq1.genes[0].to_tbl.return_value = "mockgene to tbl"
        tbl = self.seq1.to_tbl()
        expected = "1\t7\tREFERENCE\n"
        expected += "\t\t\tPBARC\t12345\n"
        expected += "mockgene to tbl"
        self.assertEquals(tbl, expected)
        
    def test_stats(self):
        self.add_mock_gene_with_1_mrna()
        self.add_mock_gene_with_2_mrnas()
        stats = self.seq1.stats()
        self.assertEquals(stats["seq_length"], 7)
        self.assertEquals(stats["num_genes"], 2)
        self.assertEquals(stats["num_mRNA"], 3)
        self.assertEquals(stats["num_CDS"], 2)



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSequence))
    return suite

if __name__ == '__main__':
    unittest.main()

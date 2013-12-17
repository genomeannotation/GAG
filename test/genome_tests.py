#!/usr/bin/env python

import unittest
from src.fasta import Fasta
from src.genome import Genome
from mock import Mock

class TestGenome(unittest.TestCase):

    def setUp(self):
        self.genome = Genome()

    def test_constructor(self):
        self.assertEqual('Genome', self.genome.__class__.__name__)

    def test_verify_file(self):
        self.assertTrue(self.genome.verify_file("sample_files/sample.sbt"))
        self.assertFalse(self.genome.verify_file("no_file_here.foo"))

    def test_remove_all_gene_segments(self):
        self.assertFalse(self.genome.gff)
        gff = Mock()
        self.genome.gff = gff
        self.assertTrue(self.genome.gff)
        self.genome.remove_all_gene_segments("BDOR_foo")
        gff.remove_all_gene_segments.assert_called_with("BDOR_foo")

    def test_remove_genes_containing_mrna_named(self):
        gff = Mock()
        self.genome.gff = gff
        self.genome.remove_genes_containing_mrna_named("foo")
        gff.remove_genes_containing_mrna_named.assert_called_with("foo")


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGenome))
    return suite

if __name__ == '__main__':
    unittest.main()

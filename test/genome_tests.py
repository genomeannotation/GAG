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

    def test_remove_mrnas_with_cds_shorter_than(self):
        gff = Mock()
        self.genome.gff = gff
        self.genome.remove_mrnas_with_cds_shorter_than(150)
        gff.remove_mrnas_with_cds_shorter_than.assert_called_with(150)

    def test_remove_first_cds_segment_if_shorter_than(self):
        gff = Mock()
        self.genome.gff = gff
        self.genome.remove_first_cds_segment_if_shorter_than(4)
        gff.remove_first_cds_segment_if_shorter_than.assert_called_with(4)

    def test_create_starts_and_stops(self):
        fasta = Mock()
        self.genome.fasta = fasta
        gene1 = Mock()
        gene2 = Mock()
        gff = Mock()
        gff.genes = [gene1, gene2]
        self.genome.gff = gff
        self.genome.create_starts_and_stops()
        gene1.create_starts_and_stops.assert_called_with(fasta)
        
    def test_obliterate_genes_related_to_mrnas(self):
        gff = Mock()
        self.genome.gff = gff
        self.genome.obliterate_genes_related_to_mrnas(['foo', 'bar'])
        gff.obliterate_genes_related_to_mrnas.assert_called_with(['foo', 'bar'])

    def setup_mocks_for_rename_maker_mrnas(self):
        self.maker_mrna = Mock()
        self.maker_mrna.name = 'maker-guy'
        self.maker_mrna.is_maker_mrna.return_value = True
        self.bdor_mrna = Mock()
        self.bdor_mrna.name = 'bdor-guy'
        self.bdor_mrna.is_maker_mrna.return_value = False
        self.gene1 = Mock()
        self.gene1.mrnas = [self.maker_mrna, self.bdor_mrna]
        self.other_mrna = Mock()
        self.other_mrna.is_maker_mrna.return_value = False
        self.other_maker = Mock()
        self.other_maker.is_maker_mrna.return_value = True
        self.gene2 = Mock()
        self.gene2.mrnas = [self.other_mrna, self.other_maker]
        self.gff = Mock()
        self.gff.genes = [self.gene1, self.gene2]
        self.genome.gff = self.gff
        self.genome.annot = Mock()

    def test_rename_maker_mrnas(self):
        self.setup_mocks_for_rename_maker_mrnas()
        self.assertEquals('maker-guy', self.maker_mrna.name)
        self.genome.rename_maker_mrnas()
        self.assertEquals('BDOR_1000000', self.maker_mrna.name)
        self.assertEquals('bdor-guy', self.bdor_mrna.name)
        self.assertEquals('BDOR_1000001', self.other_maker.name)
        self.bdor_mrna.is_maker_mrna.assert_called_with()
        self.other_mrna.is_maker_mrna.assert_called_with()

    def test_invalidate_region(self):
        gff = Mock()
        self.genome.gff = gff
        self.assertTrue(self.genome.gff)
        self.genome.invalidate_region('Scaffold_foo', 50, 100)
        gff.invalidate_region.assert_called_with('Scaffold_foo', 50, 100)

    def test_trim_region(self):
        gene1 = Mock()
        gene2 = Mock()
        gff = Mock()
        gff.genes = [gene1, gene2]
        self.genome.gff = gff
        fasta = Mock()
        self.genome.fasta = fasta
        self.genome.trim_region('seq1', 1, 3)
        fasta.trim_seq.assert_called_with('seq1', 1, 3)
        gene1.adjust_indices.assert_called_with(-3, 3) 



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGenome))
    return suite

if __name__ == '__main__':
    unittest.main()

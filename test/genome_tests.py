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

    def make_mock_genes_and_gff(self):
        self.gene1 = Mock()
        self.gene2 = Mock()
        self.gff = Mock()
        self.gff.genes = [self.gene1, self.gene2]
        self.genome.gff = self.gff

    def test_trim_region(self):
        self.make_mock_genes_and_gff()
        self.gene1.seq_name = 'seq1'
        fasta = Mock()
        self.genome.fasta = fasta
        self.genome.trim_region('seq1', 1, 3)
        fasta.trim_seq.assert_called_with('seq1', 1, 3)
        self.gene1.adjust_indices.assert_called_with(-3, 3) 

    def test_remove_seq(self):
        self.genome.fasta = Mock()
        gff = Mock()
        gff.contains_gene_on_seq.return_value = False
        self.genome.gff = gff
        self.genome.remove_seq('seq1')
        self.genome.fasta.remove_seq.assert_called_with('seq1')

    def test_get_gene_seq_info(self):
        self.make_mock_genes_and_gff()
        self.gene1.seq_name = 'seq1'
        self.gene1.indices = [1, 50]
        self.gene1.name = 'BDOR_FOO'
        self.gene2.seq_name = 'seq1'
        self.gene2.indices = [100, 150]
        self.gene2.name = 'BDOR_BAR'
        actual = self.genome.get_gene_seq_info('BDOR_FOO')
        expected = ['seq1', [1, 50]]
        self.assertEqual(actual, expected)

    def setup_mock_genes_with_name_indices_seq_name(self):
        self.make_mock_genes_and_gff()
        self.gene1.seq_name = 'seq1'
        self.gene1.indices = [1, 50]
        self.gene1.name = 'BDOR_FOO'
        self.gene2.seq_name = 'seq1'
        self.gene2.indices = [100, 150]
        self.gene2.name = 'BDOR_BAR'

    def test_check_gene_for_invalid_begin_or_end(self):
        self.setup_mock_genes_with_name_indices_seq_name()
        fasta = Mock()
        self.genome.fasta = fasta
        self.genome.fasta.how_many_Ns_forward.return_value = 5
        self.genome.fasta.how_many_Ns_backward.return_value = 0
        self.genome.check_gene_for_invalid_begin_or_end('BDOR_FOO')
        self.gff.invalidate_region.assert_called_with('seq1', 1, 5)

    def test_check_gene_for_invalid_begin_or_end_invalid_end(self):
        self.setup_mock_genes_with_name_indices_seq_name()
        fasta = Mock()
        self.genome.fasta = fasta
        self.genome.fasta.how_many_Ns_forward.return_value = 0
        self.genome.fasta.how_many_Ns_backward.return_value = 10
        self.genome.check_gene_for_invalid_begin_or_end('BDOR_BAR')
        self.gff.invalidate_region.assert_called_with('seq1', 141, 150)



        





##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGenome))
    return suite

if __name__ == '__main__':
    unittest.main()

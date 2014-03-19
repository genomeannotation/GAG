#!/usr/bin/env python

import unittest
from mock import Mock, patch, PropertyMock
import sys
import os
from src.console_controller import ConsoleController

class TestConsoleController(unittest.TestCase):

    def setUp(self):
        self.ctrlr = ConsoleController()

    def test_constructor(self):
        self.assertEqual('ConsoleController', self.ctrlr.__class__.__name__)

    def test_status(self):
        expected = "Fasta: no fasta\nGenes: no genes\n"
        self.assertEquals(expected, self.ctrlr.status())

        expected2 = "Fasta: Fasta containing 5 sequences\nGenes: 2 genes\n"
        def fastastring(self):
            return "Fasta containing 5 sequences\n"
        mock_genome = Mock()
        mock_fasta = Mock()
        mock_fasta.__str__ = fastastring
        type(mock_genome).fasta = mock_fasta
        mock_genome.genes = [Mock(), Mock()]
        self.ctrlr.genome = mock_genome
        self.assertEquals(expected2, self.ctrlr.status())

    def test_read_fasta(self):
        self.assertFalse(self.ctrlr.genome.fasta)
        self.ctrlr.read_fasta("walkthrough/gag.fasta")
        self.assertTrue(self.ctrlr.genome.fasta)

    def test_read_gff(self):
        self.assertFalse(self.ctrlr.genome.genes)
        self.ctrlr.read_gff("walkthrough/gag.gff")
        self.assertTrue(self.ctrlr.genome.genes)

    def test_ducttape(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.ducttape()
        mock_genome.rename_maker_mrnas.assert_called_with()
        mock_genome.remove_first_cds_segment_if_shorter_than.assert_called_with(3)
        mock_genome.create_starts_and_stops.assert_called_with()
        mock_genome.remove_mrnas_with_cds_shorter_than.assert_called_with(150)

    def test_create_starts_and_stops(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.create_starts_and_stops()
        mock_genome.create_starts_and_stops.assert_called_with()

    def test_subset_genome(self):
        mock_genome = Mock()
        mock_fasta = Mock()
        mock_genome.fasta = mock_fasta
        self.ctrlr.genome = mock_genome
        self.ctrlr.subset_genome('seq_foo')
        mock_genome.fasta.subset_fasta.assert_called_with(['seq_foo'])
        mock_genome.subset_genes.assert_called_with(['seq_foo'])

    def test_remove_gene(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.remove_gene("BDOR_foo")
        mock_genome.remove_all_gene_segments.assert_called_with("BDOR_foo")

    def test_remove_all_gene_segments_multiline(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.remove_gene("BDOR_foo\nBDOR_bar\nBDOR_sandwich")
        calls = mock_genome.mock_calls

        # build a list of arguments to all calls to mock_genome
        argslist = []
        for call in calls:
            argslist.append(call[1][0])

        self.assertTrue("BDOR_foo" in argslist)
        self.assertTrue("BDOR_bar" in argslist)
        self.assertTrue("BDOR_sandwich" in argslist)

    def test_rename_maker_mrnas(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.rename_maker_mrnas()
        mock_genome.rename_maker_mrnas.assert_called_with()

    def test_remove_mrnas_with_cds_shorter_than(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.remove_mrnas_with_cds_shorter_than(150)
        mock_genome.remove_mrnas_with_cds_shorter_than.assert_called_with(150)

    def test_trim_region(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.trim_region('seq1 5 8')
        mock_genome.trim_region.assert_called_with('seq1', 5, 8)

    def test_remove_seq(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.remove_seq('seq1')
        mock_genome.remove_seq.assert_called_with('seq1')

    def test_check_gene_for_invalid_begin_or_end(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.check_gene_for_invalid_begin_or_end('BDOR_FOO')
        mock_genome.check_gene_for_invalid_begin_or_end.assert_called_with('BDOR_FOO')

    def test_invalidate_region(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.invalidate_region('seq1 5 10')
        mock_genome.invalidate_region.assert_called_with('seq1', 5, 10)



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestConsoleController))
    return suite

if __name__ == '__main__':
    unittest.main()

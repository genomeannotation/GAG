#!/usr/bin/env python

import unittest
from mock import Mock, patch, PropertyMock
import sys
import os
from src.console_controller import ConsoleController, format_list_with_strings
from src.sequence import Sequence

class TestConsoleController(unittest.TestCase):

    def setUp(self):
        self.ctrlr = ConsoleController()

    def test_constructor(self):
        self.assertEqual('ConsoleController', self.ctrlr.__class__.__name__)

    def test_format_list_with_strings(self):
        mylist = ["one", "two", "three"]
        expected = "one, two, three\n"
        self.assertEquals(format_list_with_strings(mylist), expected)
    
    def test_format_list_with_strings_one_entry(self):
        self.assertEquals(format_list_with_strings(["foo"]), "foo\n")

    def test_format_list_with_strings_two_entries(self):
        self.assertEquals(format_list_with_strings(["foo", "bar"]), "foo, bar\n")

    def setup_seqs(self):
        self.ctrlr.seqs.append(Sequence("seq1", "GATTACA"))
        self.ctrlr.seqs.append(Sequence("seq2", "ATTAC"))
        self.ctrlr.seqs.append(Sequence("seq3", "ACGTACGT"))

    def setup_seqs_and_genes(self):
        self.setup_seqs()
        gene1 = Mock()
        gene1.seq_name = "seq2"
        gene1.identifier = "gene1"
        self.ctrlr.add_gene(gene1)
        gene2 = Mock()
        gene2.seq_name = "seq2"
        gene2.identifier = "gene2"
        self.ctrlr.add_gene(gene2)
        gene3 = Mock()
        gene3.seq_name = "seq3"
        gene3.identifier = "gene3"
        self.ctrlr.add_gene(gene3)

    def test_clear_seqs(self):
        self.setup_seqs()
        self.assertEquals(3, len(self.ctrlr.seqs))
        self.ctrlr.clear_seqs()
        self.assertEquals(0, len(self.ctrlr.seqs))

    def test_get_n_seq_ids(self):
        self.setup_seqs()
        expected = "First 3 seq ids are: seq1, seq2, seq3\n"
        self.assertEquals(self.ctrlr.get_n_seq_ids(3), expected)

    def test_get_n_seq_ids_when_n_is_larger_than_num_seqs(self):
        self.setup_seqs()
        expected = "First 3 seq ids are: seq1, seq2, seq3\n"
        self.assertEquals(self.ctrlr.get_n_seq_ids(8), expected)

    def test_get_n_seq_ids_when_no_seqs(self):
        expected = "No sequences currently in memory.\n"
        self.assertEquals(self.ctrlr.get_n_seq_ids(8), expected)

    def test_add_gene(self):
        self.setup_seqs()
        gene1 = Mock()
        gene1.seq_name = "seq1"
        self.assertEquals(0, len(self.ctrlr.seqs[0].genes))
        self.ctrlr.add_gene(gene1)
        self.assertEquals(1, len(self.ctrlr.seqs[0].genes))

    def test_read_fasta(self):
        self.assertFalse(self.ctrlr.seqs)
        self.ctrlr.read_fasta("walkthrough/gag.fasta")
        self.assertTrue(self.ctrlr.seqs)

    def test_read_gff(self):
        self.ctrlr.read_fasta("walkthrough/gag.fasta")
        self.assertFalse(self.ctrlr.seqs[0].genes)
        self.ctrlr.read_gff("walkthrough/gag.gff")
        self.assertTrue(self.ctrlr.seqs[0].genes)

    def test_create_starts_and_stops(self):
        pass
        #mock_genome = Mock()
        #self.ctrlr.genome = mock_genome
        #self.ctrlr.create_starts_and_stops()
        #mock_genome.create_starts_and_stops.assert_called_with()

    def test_subset_genome(self):
        self.setup_seqs()
        self.assertEquals(3, len(self.ctrlr.seqs))
        self.ctrlr.subset_genome("seq1 seq3")
        self.assertEquals(2, len(self.ctrlr.seqs))

    def test_remove_gene(self):
        self.setup_seqs_and_genes()
        self.assertEquals(2, len(self.ctrlr.seqs[1].genes))
        self.ctrlr.remove_gene("gene2")
        self.assertEquals(1, len(self.ctrlr.seqs[1].genes))

    def test_remove_all_genes_multiline(self):
        self.setup_seqs_and_genes()
        self.assertEquals(2, len(self.ctrlr.seqs[1].genes))
        self.assertEquals(1, len(self.ctrlr.seqs[2].genes))
        self.ctrlr.remove_gene("gene2\ngene1")
        self.assertEquals(0, len(self.ctrlr.seqs[1].genes))
        self.assertEquals(1, len(self.ctrlr.seqs[2].genes))

    def test_trim_region(self):
        # TODO this isn't quite right ... gene.trim needs looking at.
        self.setup_seqs_and_genes()
        self.assertEquals(8, len(self.ctrlr.seqs[2].bases))
        #self.ctrlr.trim_region("seq3 1 3")
        #self.assertEquals(5, len(self.ctrlr.seqs[2].bases))
        # TODO after add Sequence.trim()
        #self.ctrlr.seqs.append(Sequence("seq1", "GATTACA"))
        #self.ctrlr.seqs.append(Sequence("seq2", "ATTAC"))
        #self.ctrlr.seqs.append(Sequence("seq3", "ACGTACGT"))

    def test_remove_seq(self):
        self.setup_seqs()
        self.assertEquals(3, len(self.ctrlr.seqs))
        self.ctrlr.remove_seq("seq2")
        self.assertEquals(2, len(self.ctrlr.seqs))

    def test_remove_seq_force_necessary(self):
        self.setup_seqs_and_genes()
        self.assertEquals(3, len(self.ctrlr.seqs))
        self.ctrlr.remove_seq("seq2")
        self.assertEquals(2, len(self.ctrlr.seqs))

    def test_check_gene_for_invalid_begin_or_end(self):
        pass
        #mock_genome = Mock()
        #self.ctrlr.genome = mock_genome
        #self.ctrlr.check_gene_for_invalid_begin_or_end('BDOR_FOO')
        #mock_genome.check_gene_for_invalid_begin_or_end.assert_called_with('BDOR_FOO')

    def test_invalidate_region(self):
        pass
        #mock_genome = Mock()
        #self.ctrlr.genome = mock_genome
        #self.ctrlr.invalidate_region('seq1 5 10')
        #mock_genome.invalidate_region.assert_called_with('seq1', 5, 10)

    def test_barfseq_no_args(self):
        pass
        line = ""
        # Shouldn't throw error
        self.assertTrue(self.ctrlr.barf_seq(line))

    def test_barfseq(self):
        self.setup_seqs()
        result = self.ctrlr.barf_seq("seq1 1 3")
        self.assertEquals("GAT", result)


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestConsoleController))
    return suite

if __name__ == '__main__':
    unittest.main()

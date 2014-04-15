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

    def test_contains_mrna(self):
        mockseq = Mock()
        mockseq.contains_mrna.return_value = True
        self.ctrlr.seqs = [mockseq]
        self.assertTrue(self.ctrlr.contains_mrna("foo_mrna"))
        mockseq.contains_mrna.assert_called_with("foo_mrna")

    def test_contains_gene(self):
        mockseq = Mock()
        mockseq.contains_gene.return_value = True
        self.ctrlr.seqs = [mockseq]
        self.assertTrue(self.ctrlr.contains_gene("foo_gene"))
        mockseq.contains_gene.assert_called_with("foo_gene")

    def test_contains_seq(self):
        self.setup_seqs()
        self.assertTrue(self.ctrlr.contains_seq("seq1"))
        self.assertFalse(self.ctrlr.contains_seq("foo_seq"))

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

    def test_get_n_gene_ids(self):
        self.ctrlr.seqs = [Mock(), Mock(), Mock()]
        self.ctrlr.seqs[0].get_gene_ids.return_value = ["gene1a", "gene1b", "gene1c"]
        self.ctrlr.seqs[1].get_gene_ids.return_value = ["gene2a", "gene2b", "gene2c"]
        self.ctrlr.seqs[2].get_gene_ids.return_value = ["gene3a", "gene3b", "gene3c"]
        expected = "First 4 gene ids are: gene1a, gene1b, gene1c, gene2a\n"
        self.assertEquals(self.ctrlr.get_n_gene_ids(4), expected)

    def test_get_n_mrna_ids(self):
        self.ctrlr.seqs = [Mock(), Mock(), Mock()]
        self.ctrlr.seqs[0].get_mrna_ids.return_value = ["mrna1a", "mrna1b", "mrna1c"]
        self.ctrlr.seqs[1].get_mrna_ids.return_value = ["mrna2a", "mrna2b", "mrna2c"]
        self.ctrlr.seqs[2].get_mrna_ids.return_value = ["mrna3a", "mrna3b", "mrna3c"]
        expected = "First 7 mrna ids are: mrna1a, mrna1b, mrna1c, mrna2a, mrna2b, mrna2c, mrna3a\n"
        self.assertEquals(self.ctrlr.get_n_mrna_ids(7), expected)

    def test_read_fasta(self):
        self.assertFalse(self.ctrlr.seqs)
        self.ctrlr.read_fasta("walkthrough/gag.fasta")
        self.assertTrue(self.ctrlr.seqs)

    def test_read_gff(self):
        self.ctrlr.read_fasta("walkthrough/gag.fasta")
        self.assertFalse(self.ctrlr.seqs[0].genes)
        self.ctrlr.read_gff("walkthrough/gag.gff")
        self.assertTrue(self.ctrlr.seqs[0].genes)

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

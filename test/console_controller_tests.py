#!/usr/bin/env python

import unittest
from mock import Mock, patch, PropertyMock
import sys
import io
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

    def test_add_annotations_from_list(self):
        self.ctrlr.seqs = [Mock(), Mock()]
        anno_list = [["gene1", "name", "ABC123"], ["foo_mrna", "Dbxref", "PFAM:0001"]]
        self.ctrlr.add_annotations_from_list(anno_list)
        self.ctrlr.seqs[0].add_annotations_from_list.assert_called_with(anno_list)
        self.ctrlr.seqs[1].add_annotations_from_list.assert_called_with(anno_list)

    def test_remove_from_list(self):
        self.setup_seqs_and_genes()
        bad_items = ["seq2", "foo_gene", "bar_mrna"]
        self.assertEquals(3, len(self.ctrlr.seqs))
        self.ctrlr.remove_from_list(bad_items)
        # Verify 'seq2' removed
        self.assertEquals(2, len(self.ctrlr.seqs))
        self.assertEquals("seq3", self.ctrlr.seqs[1].header)

    def test_trim_from_list_beginning_of_seq(self):
        self.setup_seqs()
        trimlist = [["seq1", 1, 3]]
        self.assertEquals(7, len(self.ctrlr.seqs[0].bases))
        self.ctrlr.trim_from_list(trimlist)
        self.assertEquals(4, len(self.ctrlr.seqs[0].bases))

    def test_trim_from_list_end_of_seq(self):
        self.setup_seqs()
        trimlist = [["seq2", 3, 5]]
        self.assertEquals(5, len(self.ctrlr.seqs[1].bases))
        self.ctrlr.trim_from_list(trimlist)
        self.assertEquals(2, len(self.ctrlr.seqs[1].bases))

    def test_trim_from_list_beginning_and_end_of_seq(self):
        self.setup_seqs()
        trimlist = [["seq3", 1, 2], ["seq3", 6, 8]]
        self.assertEquals(8, len(self.ctrlr.seqs[2].bases))
        self.ctrlr.trim_from_list(trimlist)
        self.assertEquals(3, len(self.ctrlr.seqs[2].bases))
        self.assertEquals("GTA", self.ctrlr.seqs[2].bases)

    def test_read_annotation_file(self):
        annofile = io.BytesIO("foo_gene\tname\tABC123\nfoo_mrna\tDbxref\tPFAM:0001\n")
        result = self.ctrlr.read_annotation_file(annofile)
        expected = [["foo_gene", "name", "ABC123"], ["foo_mrna", "Dbxref", "PFAM:0001"]]
        self.assertEquals(expected, result)

    def test_read_bed_file(self):
        trimfile = io.BytesIO("seq3\t1\t2\nseq3\t6\t8\n")
        result = self.ctrlr.read_bed_file(trimfile)
        expected = [["seq3", 1, 2], ["seq3", 6, 8]]
        self.assertEquals(expected, result)

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
        self.ctrlr.read_fasta("walkthrough/basic/genome.fasta")
        self.assertTrue(self.ctrlr.seqs)

    def test_read_gff(self):
        self.ctrlr.read_fasta("walkthrough/basic/genome.fasta")
        self.assertFalse(self.ctrlr.seqs[0].genes)
        self.ctrlr.read_gff("walkthrough/basic/genome.gff")
        self.assertTrue(self.ctrlr.seqs[0].genes)

    def test_barfseq_no_args(self):
        pass
        line = ""
        # Shouldn't throw error
        self.assertTrue(self.ctrlr.barf_seq(line))

    def test_barfseq(self):
        self.setup_seqs()
        result = self.ctrlr.barf_seq("seq1 1 3")
        self.assertEquals("GAT", result)

    def test_can_write_to_path(self):
        self.assertFalse(self.ctrlr.can_write_to_path("src/"))
        self.assertFalse(self.ctrlr.can_write_to_path("gag.py"))
        self.assertTrue(self.ctrlr.can_write_to_path("no_such_directory/no_such_subdirectory/no_such_file.txt"))



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestConsoleController))
    return suite

if __name__ == '__main__':
    unittest.main()

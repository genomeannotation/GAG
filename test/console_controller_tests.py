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
        self.assertTrue(isinstance(self.ctrlr.seqlist, list))

    def test_status(self):
        expected = "Fasta: no fasta\nGFF: no gff\nTemplate File: no template file\n"
        expected += "Seqlist: no seqlist\nTbl2asn Executable: no tbl2asn executable\n"
        self.assertEquals(expected, self.ctrlr.status())

        expected2 = "Fasta: Fasta containing 5 sequences\nGFF: GFF containing 20 genes\n"
        expected2 += "Template File: foo.sbt\nSeqlist: ['seq1', 'seq3']\n"
        expected2 += "Tbl2asn Executable: linux.tbl2asn\n"
        def fastastring(self):
            return "Fasta containing 5 sequences\n"
        def gffstring(self):
            return "GFF containing 20 genes\n"
        mock_genome = Mock()
        mock_fasta = Mock()
        mock_fasta.__str__ = fastastring
        type(mock_genome).fasta = mock_fasta
        mock_gff = Mock()
        mock_gff.__str__ = gffstring
        type(mock_genome).gff = mock_gff
        self.ctrlr.genome = mock_genome
        self.ctrlr.template_file = "foo.sbt"
        self.ctrlr.seqlist = ['seq1', 'seq3']
        self.ctrlr.tbl2asn_executable = "linux.tbl2asn"
        self.assertEquals(expected2, self.ctrlr.status())

    def test_read_fasta(self):
        self.assertFalse(self.ctrlr.genome.fasta)
        self.ctrlr.read_fasta("demo/demo.fasta")
        self.assertTrue(self.ctrlr.genome.fasta)

    def test_add_seq(self):
        self.assertEquals(0, len(self.ctrlr.seqlist))
        self.ctrlr.add_seq('seq2')
        self.assertEquals(1, len(self.ctrlr.seqlist))

    def test_add_seq_multiline(self):
        self.assertEquals(0, len(self.ctrlr.seqlist))
        self.ctrlr.add_seq("seq1\nseq2\nseq3\nseq4\n")
        self.assertEquals(4, len(self.ctrlr.seqlist))
        self.assertEquals('seq3', self.ctrlr.seqlist[2])

    def test_clear_seqlist(self):
        self.assertEquals(0, len(self.ctrlr.seqlist))
        self.ctrlr.seqlist = ['fooseq', 'barseq']
        self.assertEquals(2, len(self.ctrlr.seqlist))
        self.ctrlr.clear_seqlist()
        self.assertEquals(0, len(self.ctrlr.seqlist))

    def test_add_template_file(self):
        self.assertFalse(self.ctrlr.template_file)
        self.ctrlr.add_template_file("demo/demo.sbt")
        self.assertTrue(self.ctrlr.template_file)

    def test_read_gff(self):
        self.assertFalse(self.ctrlr.genome.gff)
        self.ctrlr.read_gff("demo/demo.gff")
        self.assertTrue(self.ctrlr.genome.gff)

    def test_prep_tbl2asn(self):
        self.assertFalse(os.path.isdir("tbl2asn_unittest"))
        # must have a fasta, a tbl and an sbt
        self.ctrlr.read_fasta("demo/demo.fasta")
        self.ctrlr.read_gff("demo/demo.gff")
        self.ctrlr.add_template_file("demo/demo.sbt")
        self.ctrlr.prep_tbl2asn("tbl2asn_unittest")
        self.assertTrue(os.path.isdir("tbl2asn_unittest"))
        self.assertTrue(os.path.exists("tbl2asn_unittest/gag.sbt"))
        self.assertTrue(os.path.exists("tbl2asn_unittest/gag.fsa"))
        self.assertTrue(os.path.exists("tbl2asn_unittest/gag.tbl"))
        os.system('rm -r tbl2asn_unittest')

   
    def test_ready_for_tbl2asn(self):
        self.ctrlr.set_tbl2asn_executable("actual/path/goes/here")
        def no_fsa(path):
            if path == "tbl2asn_demo/gag.fsa":
                return False
            else:
                return True
        def no_tbl(path):
            if path == "tbl2asn_demo/gag.tbl":
                return False
            else:
                return True
        def no_sbt(path):
            if path == "tbl2asn_demo/gag.sbt":
                return False
            else:
                return True
        self.assertTrue(self.ctrlr.ready_for_tbl2asn('tbl2asn_demo'))
        mock1 = Mock(side_effect=no_fsa)
        with patch('os.path.exists', mock1):
            self.assertFalse(self.ctrlr.ready_for_tbl2asn('tbl2asn_demo'))
        mock2 = Mock(side_effect=no_tbl)
        with patch('os.path.exists', mock2):
            self.assertFalse(self.ctrlr.ready_for_tbl2asn('tbl2asn_demo'))
        mock3 = Mock(side_effect=no_sbt)
        with patch('os.path.exists', mock3):
            self.assertFalse(self.ctrlr.ready_for_tbl2asn('tbl2asn_demo'))
        self.assertTrue(self.ctrlr.ready_for_tbl2asn('tbl2asn_demo'))

    def test_remove_all_gene_segments(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.remove_all_gene_segments("BDOR_foo")
        mock_genome.remove_all_gene_segments.assert_called_with("BDOR_foo")

    def test_remove_all_gene_segments_multiline(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.remove_all_gene_segments("BDOR_foo\nBDOR_bar\nBDOR_sandwich")
        calls = mock_genome.mock_calls

        # build a list of arguments to all calls to mock_genome
        argslist = []
        for call in calls:
            argslist.append(call[1][0])

        self.assertTrue("BDOR_foo" in argslist)
        self.assertTrue("BDOR_bar" in argslist)
        self.assertTrue("BDOR_sandwich" in argslist)

    def test_obliterate_genes_related_to_mrnas_multiline(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.obliterate_genes_related_to_mrnas("BDOR_foo\nBDOR_bar\nBDOR_sandwich")
        mock_genome.obliterate_genes_related_to_mrnas.assert_called_with(["BDOR_foo", "BDOR_bar", "BDOR_sandwich"])

    def test_rename_maker_mrnas(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.rename_maker_mrnas()
        mock_genome.rename_maker_mrnas.assert_called_with()

    def test_verify_all_starts_and_stops(self):
        mock_genome = Mock()
        self.ctrlr.genome = mock_genome
        self.ctrlr.verify_all_starts_and_stops()
        mock_genome.verify_all_starts_and_stops.assert_called_with()




##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestConsoleController))
    return suite

if __name__ == '__main__':
    unittest.main()

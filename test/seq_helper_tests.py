#!/usr/bin/env python
# coding=utf-8

import unittest

from mock import Mock

from src.seq_helper import SeqHelper


class TestSeqHelper(unittest.TestCase):
    def setUp(self):
        self.helper = SeqHelper("nnnGATTACAnnnnnNNNNNgattacaNNN")
        self.mrna = Mock()
        self.mrna.identifier = "foo_mrna_identifier"
        # ToDo: check if parent_id and name can be assumed to exist
        self.mrna.parent_id = "foo_mrna_parent_id"
        self.mrna.name = "foo_mrna_name"

    def tearDown(self):
        del self.mrna
        del self.helper

    def test_initialize(self):
        expected = "nnnGATTACAnnnnnNNNNNgattacaNNN"
        calculated = self.helper.full_sequence
        self.assertEquals(expected, calculated)

    def test_mrna_to_fasta(self):
        self.mrna.exon = Mock()
        self.mrna.exon.indices = [[4, 10], [21, 27]]
        self.mrna.strand = '+'
        # ToDo: check with Scott if this new syntax is acceptable/desired
        # expected = ">foo_mrna_identifier\nGATTACAgattaca\n"
        expected = ('>foo_mrna_identifier ID=foo_mrna_identifier|Parent=foo_mrna_parent_id|Name=foo_mrna_name\n'
                    'GATTACAgattaca\n')
        calculated = self.helper.mrna_to_fasta(self.mrna)
        self.assertEquals(expected, calculated)

    def test_mrna_to_fasta_reverse_strand(self):
        self.mrna.exon = Mock()
        self.mrna.exon.indices = [[4, 10], [21, 27]]
        self.mrna.strand = '-'
        # ToDo: check with Scott if this new syntax is acceptable/desired
        # expected = ">foo_mrna_identifier\ntgtaatcTGTAATC\n"
        expected = ('>foo_mrna_identifier ID=foo_mrna_identifier|Parent=foo_mrna_parent_id|Name=foo_mrna_name\n'
                    'tgtaatcTGTAATC\n')
        calculated = self.helper.mrna_to_fasta(self.mrna)
        self.assertEquals(expected, calculated)

    def test_mrna_to_cds_fasta(self):
        self.mrna.cds = Mock()
        self.mrna.cds.indices = [[4, 10], [21, 27]]
        self.mrna.strand = '-'
        # ToDo: check with Scott if this new syntax is acceptable/desired
        # expected = ">foo_mrna_identifier CDS\ntgtaatcTGTAATC\n"
        expected = ('>CDS|foo_mrna_identifier ID=foo_mrna_identifier|Parent=foo_mrna_parent_id|Name=foo_mrna_name\n'
                    'tgtaatcTGTAATC\n')
        calculated = self.helper.mrna_to_cds_fasta(self.mrna)
        self.assertEquals(expected, calculated)

    def test_mrna_to_protein_fasta(self):
        self.mrna.cds = Mock()
        self.mrna.cds.indices = [[4, 10], [21, 27]]
        self.mrna.cds.get_phase.return_value = 0
        self.mrna.strand = '+'
        # ToDo: check with Scott if this new syntax is acceptable/desired
        # expected = ">foo_mrna_identifier protein\nDYRL\n"
        expected = ('>protein|foo_mrna_identifier ID=foo_mrna_identifier|Parent=foo_mrna_parent_id|Name=foo_mrna_name\n'
                    'DYRL\n')
        calculated = self.helper.mrna_to_protein_fasta(self.mrna)
        self.assertEquals(expected, calculated)

    def test_mrna_to_protein_fasta_phase1(self):
        self.mrna.cds = Mock()
        self.mrna.cds.indices = [[4, 10], [21, 27]]
        self.mrna.cds.get_phase.return_value = 1
        self.mrna.strand = '+'
        # ToDo: check with Scott if this new syntax is acceptable/desired
        # expected = ">foo_mrna_identifier protein\nITDY\n"
        expected = ('>protein|foo_mrna_identifier ID=foo_mrna_identifier|Parent=foo_mrna_parent_id|Name=foo_mrna_name\n'
                    'ITDY\n')
        calculated = self.helper.mrna_to_protein_fasta(self.mrna)
        self.assertEquals(expected, calculated)

    def test_mrna_to_protein_fasta_reverse(self):
        self.mrna.cds = Mock()
        self.mrna.cds.indices = [[21, 27], [4, 10]]
        self.mrna.cds.get_phase.return_value = 0
        self.mrna.strand = '-'
        # ToDo: check with Scott if this new syntax is acceptable/desired
        # expected = ">foo_mrna_identifier protein\nCNL*\n"
        expected = ('>protein|foo_mrna_identifier ID=foo_mrna_identifier|Parent=foo_mrna_parent_id|Name=foo_mrna_name\n'
                    'CNL*\n')
        calculated = self.helper.mrna_to_protein_fasta(self.mrna)
        self.assertEquals(expected, calculated)

    def test_mrna_to_protein_fasta_reverse_phase1(self):
        self.mrna.cds = Mock()
        self.mrna.cds.indices = [[21, 27], [4, 10]]
        self.mrna.cds.get_phase.return_value = 1
        self.mrna.strand = '-'
        # ToDo: check with Scott if this new syntax is acceptable/desired
        # expected = ">foo_mrna_identifier protein\nVICN\n"
        expected = ('>protein|foo_mrna_identifier ID=foo_mrna_identifier|Parent=foo_mrna_parent_id|Name=foo_mrna_name\n'
                    'VICN\n')
        calculated = self.helper.mrna_to_protein_fasta(self.mrna)
        self.assertEquals(expected, calculated)

    def test_get_sequence_from_indices(self):
        self.mrna.cds = Mock()
        self.mrna.cds.indices = [[21, 27], [4, 10]]
        self.mrna.strand = '+'
        expected = 'gattacaGATTACA'
        calculated = self.helper.get_sequence_from_indices(self.mrna.strand, self.mrna.cds.indices)
        self.assertEquals(expected, calculated)

    def test_get_sequence_from_indices_reverse(self):
        self.mrna.cds = Mock()
        self.mrna.cds.indices = [[21, 27], [4, 10]]
        self.mrna.strand = '-'
        expected = 'TGTAATCtgtaatc'
        calculated = self.helper.get_sequence_from_indices(self.mrna.strand, self.mrna.cds.indices)
        self.assertEquals(expected, calculated)

    def test_mrna_contains_internal_stop(self):
        helper = SeqHelper("gattacaTAGgattaca")  # TAG = stop codon
        self.mrna.strand = '+'
        self.mrna.cds = Mock()
        self.mrna.cds.indices = [[2, 4], [8, 14]]
        # verify that the internal stop is in the sequence ...
        expected = "attTAGgatt"
        calculated = helper.get_sequence_from_indices('+', self.mrna.cds.indices)
        self.assertEquals(expected, calculated)
        # make sure the method catches it
        test = helper.mrna_contains_internal_stop(self.mrna)
        self.assertTrue(test)


def suite():
    _suite = unittest.TestSuite()
    _suite.addTest(unittest.makeSuite(TestSeqHelper))
    return _suite


if __name__ == '__main__':
    unittest.main()

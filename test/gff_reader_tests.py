#!/usr/bin/env python

import unittest
import io
from mock import Mock, patch, PropertyMock
from src.gff_reader import *
import sys
import os
import csv

class TestGFFReader(unittest.TestCase):

    def setUp(self):
        self.reader = GFFReader()

    def test_validate_line_not_enough_fields(self):
        badline = "scaffold00080\tmaker\tgene\t106151\t109853\t+\t.\tID=BDOR_007864\n"
        self.assertFalse(self.reader.validate_line(badline))

    def test_validate_line_no_id(self):
        badline = "scaffold00080\tmaker\tgene\t106151\t109853\t.\t+\t.\tName=BDOR_007864\n"
        self.assertFalse(self.reader.validate_line(badline))

    def test_validate_line(self):
        goodline = "scaffold00080\tmaker\tgene\t106151\t109853\t.\t+\t.\tID=BDOR_007864\n"
        self.assertTrue(self.reader.validate_line(goodline))

    def test_line_type_gene(self):
        line = "scaffold00080\tmaker\tgene\t106151\t109853\t.\t+\t.\tID=BDOR_007864\n".split('\t')
        self.assertEqual('gene', self.reader.line_type(line))

    def test_line_type_mrna(self):
        line = "scaffold00080\tmaker\tmRNA\t106151\t109853\t.\t+\t.\tID=BDOR_007864-RA;Parent=BDOR_007864\n".split('\t')
        self.assertEqual('mRNA', self.reader.line_type(line))

    def test_line_type_exon(self):
        line = "scaffold00080\tmaker\texon\t106151\t106451\t0.9\t+\t.\tID=BDOR_007864-RA:exon:0;Parent=BDOR_007864-RA\n".split('\t')
        self.assertEqual('exon', self.reader.line_type(line))

    def test_line_type_cds(self):
        line = "scaffold00080\tmaker\tCDS\t106151\t106451\t.\t+\t0\tID=BDOR_007864-RA:cds:0;Parent=BDOR_007864-RA\n".split('\t')
        self.assertEqual('CDS', self.reader.line_type(line))

    def test_line_type_start_codon(self):
        line = "scaffold00080\tmaker\tstart_codon\t106151\t106153\t.\t+\t.\tID=BDOR_007864-RA:start1;Parent=BDOR_007864-RA\n".split('\t')
        self.assertEqual('start_codon', self.reader.line_type(line))

    def test_line_type_stop_codon(self):
        line = "scaffold00080\tmaker\tstop_codon\t109851\t109853\t.\t+\t.\tID=BDOR_007864-RA:stop2;Parent=BDOR_007864-RA\n".split('\t')
        self.assertEqual('stop_codon', self.reader.line_type(line))

    def test_parse_attributes(self):
        attr = "ID=BDOR_007864-RA:stop2;Parent=BDOR_007864-RA\n"
        parsed = self.reader.parse_attributes(attr)
        self.assertEqual('BDOR_007864-RA:stop2', parsed['identifier'])
        self.assertEqual('BDOR_007864-RA', parsed['parent_id'])

    def test_extract_cds_args(self):
        line = "scaffold00080\tmaker\tCDS\t106151\t106451\t.\t+\t0\tID=BDOR_007864-RA:cds:0;Parent=BDOR_007864-RA\n".split('\t')
        args = self.reader.extract_cds_args(line)
        expected = {'indices': [106151, 106451], 'phase': 0, 'identifier': 'BDOR_007864-RA:cds:0', 'parent_id': 'BDOR_007864-RA'}
        self.assertEqual(expected, args)

    def test_extract_exon_args(self):
        line = "scaffold00080\tmaker\texon\t106151\t106451\t0.9\t+\t.\tID=BDOR_007864-RA:exon:0;Parent=BDOR_007864-RA\n".split('\t')
        expected = {'indices': [106151, 106451], 'score': 0.9, 'identifier': 'BDOR_007864-RA:exon:0', 'parent_id': 'BDOR_007864-RA'}
        args = self.reader.extract_exon_args(line)
        self.assertEqual(expected, args)

    def test_extract_mrna_args(self):
        line = "scaffold00080\tmaker\tmRNA\t106151\t109853\t.\t+\t.\tID=BDOR_007864-RA;Parent=BDOR_007864\n".split('\t')
        expected = {'indices': [106151, 109853], 'identifier': 'BDOR_007864-RA', 'parent_id': 'BDOR_007864'}
        args = self.reader.extract_mrna_args(line)
        self.assertEqual(expected, args)

    def test_extract_gene_args(self):
        line = "scaffold00080\tmaker\tgene\t106151\t109853\t.\t+\t.\tID=BDOR_007864\n".split('\t')
        expected = {'seq_name': 'scaffold00080', 'source': 'maker', 'indices': [106151, 109853],\
                    'strand': '+', 'identifier': 'BDOR_007864'}
        args = self.reader.extract_gene_args(line)
        self.assertEqual(expected, args)

    def test_process_line(self):
        pass

    def test_process_all_lines(self):
        pass
        
    def get_sample_text(self):
        sample_text = "scaffold00080\tmaker\tgene\t106151\t109853\t.\t+\t.\tID=BDOR_007864\n"
        sample_text += "scaffold00080\tmaker\tmRNA\t106151\t109853\t.\t+\t.\tID=BDOR_007864-RA;Parent=BDOR_007864\n"
        sample_text += "scaffold00080\tmaker\texon\t106151\t106451\t0\t9	+\t.\tID=BDOR_007864-RA:exon:0;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\texon\t106509\t106749\t0\t9	+\t.\tID=BDOR_007864-RA:exon:1;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tCDS\t106151\t106451\t.\t+\t0	I\t=BDOR_007864-RA:cds:0;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tCDS\t106509\t106749\t.\t+\t2	I\t=BDOR_007864-RA:cds:1;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tstart_codon\t106151\t106153\t.\t+\t.\tID=BDOR_007864-RA:start1;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tstop_codon\t109851\t109853\t.\t+\t.\tID=BDOR_007864-RA:stop2;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tgene\t145206\t183302\t.\t+\t.\tID=BDOR_007866\n"
        sample_text += "scaffold00080\tmaker\tmRNA\t145206\t183302\t.\t+\t.\tID=BDOR_007866-RB;Parent=BDOR_007866\n"
        sample_text += "scaffold00080\tmaker\texon\t145206\t145282\t0\t065333	+\t.\tID=BDOR_007866-RB:exon:5;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\texon\t145607\t145865\t47\t919	+\t.\tID=BDOR_007866-RB:exon:6;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\texon\t145928\t146176\t67\t378	+\t.\tID=BDOR_007866-RB:exon:7;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tfive_prime_UTR\t145206\t145282\t.\t+\t.\tID=BDOR_007866-RB:UTR1;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tfive_prime_UTR\t145607\t145865\t.\t+\t.\tID=BDOR_007866-RB:UTR2;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tfive_prime_UTR\t145928\t146176\t.\t+\t.\tID=BDOR_007866-RB:UTR3;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tfive_prime_UTR\t154498\t154575\t.\t+\t.\tID=BDOR_007866-RB:UTR4;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tCDS\t154576\t154620\t.\t+\t0	I\t=BDOR_007866-RB:cds:5;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tCDS\t179210\t179419\t.\t+\t0	I\t=BDOR_007866-RB:cds:6;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tCDS\t179489\t179691\t.\t+\t0	I\t=BDOR_007866-RB:cds:7;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tthree_prime_UTR\t183025\t183302\t.\t+\t.\tID=BDOR_007866-RB:UTR5;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tstop_codon\t183022\t183024\t.\t+\t.\tID=BDOR_007866-RB:stop2;Parent=BDOR_007866-RB\n"
        return sample_text

        
##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGFFReader))
    return suite

if __name__ == '__main__':
    unittest.main()

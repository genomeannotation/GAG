#!/usr/bin/env python

import unittest
import io
from mock import Mock, patch, PropertyMock
from src.gff import GFF
import sys
import os
import csv

class TestGFF(unittest.TestCase):

    def setUp(self):
        self.test_gff1 = GFF()
        self.test_gff2 = GFF()
        self.test_gene_line1 = ['sctg_0080_0020', 'maker', 'gene', '3734', '7436', '.', '+', '.', 'ID=1;Name=BDOR_007864']
        self.test_gene_line2 = ['sctg_0080_0020', 'maker', 'gene', '460', '12713', '.', '-', '.', 'ID=172.1;Name=BDOR_007863.1']
        self.test_line2 = ['sctg_0080_0020', 'maker', 'mRNA', '3734', '7436', '.', '+', '.', 'ID=2;Name=BDOR_007864-RA;Parent=1']
        self.test_cds_line1 = ['sctg_0080_0020', 'maker', 'CDS', '10247', '10625', '.', '-', '1', 'ID=229.1;Name=BDOR_007863.1-RA:cds:44;Parent=173.1']
        self.test_cds_line2 = ['sctg_0080_0020', 'maker', 'CDS', '9089', '10086', '.', '-', '2', 'ID=230.1;Name=BDOR_007863.1-RA:cds:45;Parent=173.1']
        self.test_exon_line1 = ['sctg_0080_0020', 'maker', 'exon', '10247', '10625', '100.148', '-', '.', 'ID=199.1;Name=BDOR_007863.1-RA:exon:48;Parent=173.1']
        self.test_exon_line2 = ['sctg_0080_0020', 'maker', 'exon', '9089', '10086', '270.928', '-', '.', 'ID=200.1;Name=BDOR_007863.1-RA:exon:49;Parent=173.1']
        self.test_feature_line = ['sctg_0080_0020', 'maker', 'five_prime_UTR', '460', '627', '.', '-', '.', 'ID=234.1;Name=BDOR_007863.1-RA:UTR1;Parent=173.1']
        self.test_mrna_line1 = ['sctg_0080_0020', 'maker', 'mRNA', '460', '12713', '.', '-', '.', 'ID=173.1;Name=BDOR_007863.1-RA;Parent=172.1']

    def test_constructor(self):
        self.assertEqual('GFF', self.test_gff1.__class__.__name__)
        self.assertEqual(0, len(self.test_gff1.genes))
        self.assertFalse(self.test_gff1.current_gene)
        self.assertFalse(self.test_gff1.current_mrna)
        self.assertFalse(self.test_gff1.current_exon)
        self.assertFalse(self.test_gff1.current_cds)

    def test_str(self):
        expected = "GFF containing 0 genes\n"
        self.assertEquals(expected, str(self.test_gff1))

    def test_validate_line(self):
        bad_line1 = ['sctg_x', 'maker']
        self.assertFalse(self.test_gff1.validate_line(bad_line1))
        self.assertTrue(self.test_gff1.validate_line(self.test_gene_line1))

    def test_line_type(self):
        self.assertEqual('gene', self.test_gff1.line_type(self.test_gene_line1))
        self.assertEqual('mRNA', self.test_gff1.line_type(self.test_line2))

    def test_parse_attributes(self):
        expected1 = {'identifier': '2', 'parent_id': '1'}
        string1 = 'ID=2;Name=BDOR_007864-RA;Parent=1'
        self.assertEqual(expected1, self.test_gff1.parse_attributes(string1))
        expected2 = {'identifier': '1'}
        string2 = 'ID=1;Name=BDOR_007864'
        self.assertEqual(expected2, self.test_gff1.parse_attributes(string2))
        bad_string1 = 'not semicolon delimited...'
        sys.stderr.write("test should generate an error message:\n")
        self.assertFalse(self.test_gff1.parse_attributes(bad_string1))
        bad_string2 = 'semicolon delimited; no equals signs though'
        sys.stderr.write("test should generate an error message:\n")
        self.assertFalse(self.test_gff1.parse_attributes(bad_string2))

    def test_extract_cds_args(self):
        expected = {'identifier': '229.1', 'indices': [10247, 10625], 'phase': 1, 'parent_id': '173.1'}
        actual = self.test_gff1.extract_cds_args(self.test_cds_line1)
        self.assertEqual(expected, actual)

    def test_extract_exon_args(self):
        expected = {'identifier': '199.1', 'indices': [10247, 10625], 'score': '100.148', 'parent_id': '173.1'}
        actual = self.test_gff1.extract_exon_args(self.test_exon_line1)
        self.assertEqual(expected, actual)

    def test_extract_mrna_args(self):
        expected = {'identifier': '173.1', 'indices': [460, 12713], 'parent_id': '172.1'}
        actual = self.test_gff1.extract_mrna_args(self.test_mrna_line1)
        self.assertEqual(expected, actual)

    def test_extract_gene_args(self):
        expected = {'seq_name': 'sctg_0080_0020', 'source': 'maker', 'indices': [3734, 7436], 'strand': '+', 'identifier': '1'}
        actual = self.test_gff1.extract_gene_args(self.test_gene_line1)
        self.assertEqual(expected, actual)

    def test_remove_first_cds_segment_if_shorter_than(self):
        gff = GFF()
        gene1 = Mock()
        gff.genes = [gene1]
        gff.remove_first_cds_segment_if_shorter_than(4)
        gene1.remove_first_cds_segment_if_shorter_than.assert_called_with(4)

    def test_remove_mrnas_with_cds_shorter_than(self):
        gff = GFF()
        gene1 = Mock()
        gene1.mrnas = None
        gff.genes = [gene1]
        gff.remove_mrnas_with_cds_shorter_than(150)
        gene1.remove_mrnas_with_cds_shorter_than.assert_called_with(150)
        self.assertEquals(0, len(gff.genes))

    def test_process_line(self):
        with patch.object(self.test_gff1, 'process_gene_line') as mock:
            self.test_gff1.process_line(self.test_gene_line1)
        mock.assert_called_with(self.test_gene_line1)
        with patch.object(self.test_gff1, 'process_mrna_line') as mock:
            self.test_gff1.process_line(self.test_mrna_line1)
        mock.assert_called_with(self.test_mrna_line1)
        with patch.object(self.test_gff1, 'process_cds_line') as mock:
            self.test_gff1.process_line(self.test_cds_line1)
        mock.assert_called_with(self.test_cds_line1)
        with patch.object(self.test_gff1, 'process_exon_line') as mock:
            self.test_gff1.process_line(self.test_exon_line1)
        mock.assert_called_with(self.test_exon_line1)
        with patch.object(self.test_gff1, 'process_other_feature_line') as mock:
            self.test_gff1.process_line(self.test_feature_line)
        mock.assert_called_with(self.test_feature_line)

    def test_process_all_lines(self):
        # This rather lengthy test walks through reading several
        # lines -- gene, gene, mrna, mrna, cds, etc.
        # Considered breaking it up into several tests
        # but that would have involved hecka mocking.
        # Hopefully what follows is readable :)

        # Start off with no genes and no current_gene
        self.assertEqual(0, len(self.test_gff2.genes))
        self.assertFalse(self.test_gff2.current_gene)
        # if no current gene, simply create new gene
        self.test_gff2.process_gene_line(self.test_gene_line1)   
        self.assertTrue(self.test_gff2.current_gene)
        # if gene already present, should wrap up current
        # gene, add it to genes[] and instantiate a new one
        self.assertEqual(0, len(self.test_gff2.genes))
        self.test_gff2.process_gene_line(self.test_gene_line2)
        self.assertEqual(1, len(self.test_gff2.genes))
        self.assertEqual('172.1', self.test_gff2.current_gene.identifier)

        # test process_mrna_line
        self.assertFalse(self.test_gff2.current_mrna)
        # if no current mrna, should create a new one
        self.test_gff2.process_mrna_line(self.test_mrna_line1) 
        self.assertTrue(self.test_gff2.current_mrna)
        # if mrna already present, should wrap up current
        # mrna, add it to current_gene, instantiate new mrna
        self.assertEqual(0, len(self.test_gff2.current_gene.mrnas))
        self.assertEqual('173.1', self.test_gff2.current_mrna.identifier)
        self.test_mrna_line2 = ['sctg_0080_0020', 'maker', 'mRNA', '460', '12713', '.', '-', '.', 'ID=foo;Name=BDOR_007863.1-RA;Parent=172.1']
        self.test_gff2.process_mrna_line(self.test_mrna_line2)
        self.assertEqual(1, len(self.test_gff2.current_gene.mrnas))
        self.assertEqual('foo', self.test_gff2.current_mrna.identifier)

        # test process_cds_line
        self.assertFalse(self.test_gff2.current_cds)
        # if no current cds should instantiate one
        self.test_gff2.process_cds_line(self.test_cds_line1)
        self.assertTrue(self.test_gff2.current_cds)
        self.assertEqual(['229.1'], self.test_gff2.current_cds.identifier)
        # if cds already present, should add info to current_cds
        self.test_gff2.process_cds_line(self.test_cds_line2)
        self.assertEqual(['229.1', '230.1'], self.test_gff2.current_cds.identifier)

        # test process_exon_line
        self.assertFalse(self.test_gff2.current_exon)
        # if no current exon should instantiate one
        self.test_gff2.process_exon_line(self.test_exon_line1)
        self.assertTrue(self.test_gff2.current_exon)
        self.assertEqual(1, len(self.test_gff2.current_exon.indices))
        # if exon already present, add info to current_exon
        self.test_gff2.process_exon_line(self.test_exon_line2)
        self.assertEqual(2, len(self.test_gff2.current_exon.indices))

        # test process_other_feature_line
        self.assertEqual(0, len(self.test_gff2.current_mrna.other_features))
        # should do NOTHING when another feature is encountered :)
        self.test_gff2.process_other_feature_line(self.test_feature_line)
        self.assertEqual(0, len(self.test_gff2.current_mrna.other_features))
        

    def test_gff_file_stuff(self):
        sample_text = "sctg_0080_0020\tmaker\tgene\t3734\t7436\t.\t+\t.\tID=1;Name=BDOR_007864\n"
        sample_text += "sctg_0080_0020\tmaker\tmRNA\t3734\t7436\t.\t+\t.\tID=2;Name=BDOR_007864-RA;Parent=1\n"
        sample_text += "sctg_0080_0020\tmaker\texon\t3734\t4034\t0.9\t+\t.\tID=3;Name=BDOR_007864-RA:exon:0;Parent=2\n"
        sample_text += "sctg_0080_0020\tmaker\texon\t4092\t4332\t0.9\t+\t.\tID=4;Name=BDOR_007864-RA:exon:1;Parent=2\n"
        sample_text += "sctg_0080_0020\tmaker\texon\t4399\t5185\t0.9\t+\t.\tID=5;Name=BDOR_007864-RA:exon:2;Parent=2\n"
        sample_text += "sctg_0080_0020\tmaker\texon\t5249\t6565\t0.9\t+\t.\tID=6;Name=BDOR_007864-RA:exon:3;Parent=2\n"
        sample_text += "sctg_0080_0020\tmaker\texon\t6630\t7436\t0.9\t+\t.\tID=7;Name=BDOR_007864-RA:exon:4;Parent=2\n"
        sample_text += "sctg_0080_0020\tmaker\tCDS\t3734\t4034\t.\t+\t0\tID=8;Name=BDOR_007864-RA:cds:0;Parent=2\n"
        sample_text += "sctg_0080_0020\tmaker\tCDS\t4092\t4332\t.\t+\t2\tID=9;Name=BDOR_007864-RA:cds:1;Parent=2\n"
        sample_text += "sctg_0080_0020\tmaker\tCDS\t4399\t5185\t.\t+\t1\tID=10;Name=BDOR_007864-RA:cds:2;Parent=2\n"
        sample_text += "sctg_0080_0020\tmaker\tCDS\t5249\t6565\t.\t+\t0\tID=11;Name=BDOR_007864-RA:cds:3;Parent=2\n"
        sample_text += "sctg_0080_0020\tmaker\tCDS\t6630\t7436\t.\t+\t0\tID=12;Name=BDOR_007864-RA:cds:4;Parent=2\n"
        sample_text += "sctg_0080_0020\tmaker\tstart_codon\t3734\t3736\t.\t+\t.\tID=13;Name=BDOR_007864-RA:start1;Parent=2\n"
        sample_text += "sctg_0080_0020\tmaker\tstop_codon\t7434\t7436\t.\t+\t.\tID=14;Name=BDOR_007864-RA:stop2;Parent=2\n"
        sample_text += "sctg_0080_0026\tmaker\tgene\t3306\t4514\t.\t+\t.\tID=15.1;Name=BDOR_007866.1\n"
        sample_text += "sctg_0080_0026\tmaker\tmRNA\t3306\t4514\t.\t+\t.\tID=16.1;Name=BDOR_007866.1-RB;Parent=15.1\n"
        sample_text += "sctg_0080_0026\tmaker\texon\t3306\t3382\t0.065333\t+\t.\tID=17.1;Name=BDOR_007866.1-RB:exon:5;Parent=16.1\n"
        sample_text += "sctg_0080_0026\tmaker\texon\t3707\t3965\t47.919\t+\t.\tID=18.1;Name=BDOR_007866.1-RB:exon:6;Parent=16.1\n"
        sample_text += "sctg_0080_0026\tmaker\texon\t4028\t4276\t67.378\t+\t.\tID=19.1;Name=BDOR_007866.1-RB:exon:7;Parent=16.1\n"
        sample_text += "sctg_0080_0026\tmaker\tfive_prime_UTR\t3306\t3382\t.\t+\t.\tID=28.1;Name=BDOR_007866.1-RB:UTR1;Parent=16.1\n"
        sample_text += "sctg_0080_0026\tmaker\tfive_prime_UTR\t3707\t3965\t.\t+\t.\tID=29.1;Name=BDOR_007866.1-RB:UTR2;Parent=16.1\n"
        sample_text += "sctg_0080_0026\tmaker\tfive_prime_UTR\t4028\t4276\t.\t+\t.\tID=30.1;Name=BDOR_007866.1-RB:UTR3;Parent=16.1\n"
        sample_text += "sctg_0080_0026\tmaker\tmRNA\t3306\t4514\t.\t+\t.\tID=43.1;Name=BDOR_007866.1-RC;Parent=15.1\n"
        sample_text += "sctg_0080_0026\tmaker\texon\t3306\t3382\t0.065333\t+\t.\tID=44.1;Name=BDOR_007866.1-RB:exon:5;Parent=43.1\n"
        sample_text += "sctg_0080_0026\tmaker\texon\t3707\t3965\t47.919\t+\t.\tID=45.1;Name=BDOR_007866.1-RB:exon:6;Parent=43.1\n"
        sample_text += "sctg_0080_0026\tmaker\texon\t4028\t4276\t67.378\t+\t.\tID=46.1;Name=BDOR_007866.1-RB:exon:7;Parent=43.1\n"
        sample_text += "sctg_0080_0026\tmaker\tfive_prime_UTR\t3306\t3382\t.\t+\t.\tID=54.1;Name=BDOR_007866.1-RC:UTR1;Parent=43.1\n"
        sample_text += "sctg_0080_0026\tmaker\tfive_prime_UTR\t3707\t3809\t.\t+\t.\tID=55.1;Name=BDOR_007866.1-RC:UTR2;Parent=43.1\n"
        sample_text += "sctg_0080_0026\tmaker\tCDS\t3810\t3965\t.\t+\t0\tID=56.1;Name=BDOR_007866.1-RC:cds:13;Parent=43.1\n"
        sample_text += "sctg_0080_0026\tmaker\tCDS\t4028\t4276\t.\t+\t0\tID=57.1;Name=BDOR_007866.1-RC:cds:14;Parent=43.1\n"
        sample_text += "sctg_0080_0026\tmaker\tstart_codon\t3810\t3812\t.\t+\t.\tID=66.1;Name=BDOR_007866.1-RC:start1;Parent=43.1\n"
        sample_text += "sctg_0080_0026\tmaker\tmRNA\t3306\t4514\t.\t+\t.\tID=158.1;Name=BDOR_007866.1-RA;Parent=15.1\n"
        sample_text += "sctg_0080_0026\tmaker\texon\t3306\t3382\t0.065333\t+\t.\tID=159.1;Name=BDOR_007866.1-RB:exon:5;Parent=158.1\n"
        sample_text += "sctg_0080_0026\tmaker\texon\t3707\t3965\t47.919\t+\t.\tID=160.1;Name=BDOR_007866.1-RB:exon:6;Parent=158.1\n"
        sample_text += "sctg_0080_0026\tmaker\texon\t4028\t4287\t65.380\t+\t.\tID=161.1;Name=BDOR_007866.1-RA:exon:20;Parent=158.1\n"
        sample_text += "sctg_0080_0026\tmaker\tfive_prime_UTR\t3306\t3382\t.\t+\t.\tID=164.1;Name=BDOR_007866.1-RA:UTR1;Parent=158.1\n"
        sample_text += "sctg_0080_0026\tmaker\tfive_prime_UTR\t3707\t3809\t.\t+\t.\tID=165.1;Name=BDOR_007866.1-RA:UTR2;Parent=158.1\n"
        sample_text += "sctg_0080_0026\tmaker\tCDS\t3810\t3965\t.\t+\t0\tID=166.1;Name=BDOR_007866.1-RC:cds:13;Parent=158.1\n"
        sample_text += "sctg_0080_0026\tmaker\tCDS\t4028\t4287\t.\t+\t0\tID=167.1;Name=BDOR_007866.1-RA:cds:16;Parent=158.1\n"
        sample_text += "sctg_0080_0026\tmaker\tstart_codon\t3810\t3812\t.\t+\t.\tID=170.1;Name=BDOR_007866.1-RA:start1;Parent=158.1\n"
        sample_text += "sctg_0080_0027\tmaker\tgene\t1\t3870\t.\t+\t.\tID=15.2;Name=BDOR_007866.2\n"
        sample_text += "sctg_0080_0027\tmaker\tmRNA\t1\t3870\t.\t+\t.\tID=16.2;Name=BDOR_007866.2-RB;Parent=15.2\n"
        sample_text += "sctg_0080_0027\tmaker\tmRNA\t1\t3870\t.\t+\t.\tID=43.2;Name=BDOR_007866.2-RC;Parent=15.2\n"
        sample_text += "sctg_0080_0027\tmaker\tmRNA\t1584\t3870\t.\t+\t.\tID=68.2;Name=BDOR_007866.2-RD;Parent=15.2\n"
        sample_text += "sctg_0080_0027\tmaker\texon\t1584\t2084\t0.372975\t+\t.\tID=77.2;Name=BDOR_007866.2-RD:exon:16;Parent=68.2\n"
        sample_text += "sctg_0080_0027\tmaker\tfive_prime_UTR\t1584\t2084\t.\t+\t.\tID=78.2;Name=BDOR_007866.2-RD:UTR1;Parent=68.2\n"
        sample_text += "sctg_0080_0027\tmaker\tmRNA\t2224\t3870\t.\t+\t.\tID=91.2;Name=BDOR_007866.2-RF;Parent=15.2\n"
        sample_text += "sctg_0080_0027\tmaker\texon\t2224\t2402\t0.372975\t+\t.\tID=100.2;Name=BDOR_007866.2-RF:exon:17;Parent=91.2\n"
        sample_text += "sctg_0080_0027\tmaker\tfive_prime_UTR\t2224\t2402\t.\t+\t.\tID=101.2;Name=BDOR_007866.2-RF:UTR1;Parent=91.2\n"
        sample_text += "sctg_0080_0027\tmaker\tmRNA\t2063\t3870\t.\t+\t.\tID=114.2;Name=BDOR_007866.2-RE;Parent=15.2\n"
        sample_text += "sctg_0080_0027\tmaker\texon\t2063\t2135\t0.372975\t+\t.\tID=123.2;Name=BDOR_007866.2-RE:exon:18;Parent=114.2\n"
        sample_text += "sctg_0080_0027\tmaker\tfive_prime_UTR\t2063\t2135\t.\t+\t.\tID=124.2;Name=BDOR_007866.2-RE:UTR1;Parent=114.2\n"
        sample_text += "sctg_0080_0027\tmaker\tmRNA\t1\t3870\t.\t+\t.\tID=158.2;Name=BDOR_007866.2-RA;Parent=15.2')"
        sample = io.BytesIO(sample_text)
        reader = csv.reader(sample, delimiter='\t')

        gff = GFF()  
        gff.read_file(reader)
        self.assertEqual(3, len(gff.genes))
        self.assertEqual(3703, gff.genes[0].length())
        self.assertEqual('2', gff.genes[0].mrnas[0].identifier)
        self.assertEqual(1, len(gff.genes[0].mrnas))
        self.assertEqual('43.1', gff.genes[1].mrnas[1].identifier)
        self.assertEqual('15.2', gff.genes[2].identifier)
        self.assertEqual('68.2', gff.genes[2].mrnas[2].identifier)
        self.assertEqual('77.2', gff.genes[2].mrnas[2].exon.identifier[0])
        self.assertEqual(2084, gff.genes[2].mrnas[2].exon.indices[0][1])

        try:
            os.system('mkdir TEST_RUN')
        except:
            pass

        with open('TEST_RUN/foo', 'wb') as outfile:
            for gene in gff.genes:  
                outfile.write(gene.to_gff())

        os.system('rm -r TEST_RUN')

    def test_subset_gff(self):
        gene1 = Mock()
        seq1 = PropertyMock(return_value = 'sctg_1')
        type(gene1).seq_name = seq1
        gene2 = Mock()
        seq2 = PropertyMock(return_value = 'sctg_2')
        type(gene2).seq_name = seq2
        gene3 = Mock()
        seq3 = PropertyMock(return_value = 'sctg_3')
        type(gene3).seq_name = seq3
        self.test_gff1.genes.extend([gene1, gene2, gene3])
        self.assertEquals(3, len(self.test_gff1.genes))
        self.test_gff1.subset_gff(['sctg_1', 'sctg_3'])
        self.assertEquals(2, len(self.test_gff1.genes))
        self.assertEquals('sctg_1', self.test_gff1.genes[0].seq_name)
        self.assertEquals('sctg_3', self.test_gff1.genes[1].seq_name)
        
    def test_remove_empty_genes(self):
        nonempty_gene = Mock()
        nonempty_gene.is_empty.return_value = False
        empty_gene = Mock()
        empty_gene.is_empty.return_value = True
        self.test_gff1.genes.append(nonempty_gene)
        self.test_gff1.genes.append(empty_gene)
        self.assertEquals(2, len(self.test_gff1.genes))
        self.test_gff1.remove_empty_genes()
        self.assertEquals(1, len(self.test_gff1.genes))

    def test_remove_gene(self):
        gene1 = Mock()
        gene2 = Mock()
        gene3 = Mock()
        gene1.identifier = "BDOR_00001"
        gene2.identifier = "BDOR_00002.1"
        gene3.identifier = "BDOR_00002.2"
        self.assertEquals(0, len(self.test_gff1.genes))
        self.test_gff1.genes.extend([gene1, gene2, gene3])
        self.assertEquals(3, len(self.test_gff1.genes))
        self.assertEquals("BDOR_00001", self.test_gff1.genes[0].identifier)
        self.test_gff1.remove_all_gene_segments("BDOR_00002")
        self.assertEquals(1, len(self.test_gff1.genes))
        self.assertEquals("BDOR_00001", self.test_gff1.genes[0].identifier)

    def test_remove_all_gene_segments_handles_empty_string(self):
        gene1 = Mock()
        gene1.identifier = "BDOR_foo"
        self.assertEquals(0, len(self.test_gff1.genes))
        self.test_gff1.genes.append(gene1)
        self.assertEquals(1, len(self.test_gff1.genes))
        self.test_gff1.remove_all_gene_segments("")
        self.assertEquals(1, len(self.test_gff1.genes))

    def test_remove_genes_by_prefixes(self):
        prefixes = ["BDOR_00001", "BDOR_00002"]
        gene1 = Mock()
        gene1.identifier = "BDOR_00001.7"
        gene2 = Mock()
        gene2.identifier = "BDOR_00002.2"
        gene3 = Mock()
        gene3.identifier = "BDOR_00008.7"
        self.test_gff1.genes.extend([gene1, gene2, gene3])
        self.assertEquals(3, len(self.test_gff1.genes))
        self.test_gff1.remove_genes_by_prefixes(prefixes)
        self.assertEquals(1, len(self.test_gff1.genes))
        self.assertEquals("BDOR_00008.7", self.test_gff1.genes[0].identifier)

    def test_invalidate_region(self):
        gene1 = Mock()
        gene1.seq_name = 'Scaffold_foo'
        gene2 = Mock()
        gene2.seq_name = 'Scaffold_dog'
        self.test_gff1.genes.extend([gene1, gene2])
        self.test_gff1.invalidate_region('Scaffold_foo', 50, 100)
        gene1.invalidate_region.assert_called_with(50, 100)
        gene1.invalidate_region.assert_called_with(50, 100)
        assert not gene2.invalidate_region.called
        assert not gene2.invalidate_region.called
    
    def test_contains_gene_on_seq(self):
        gene1 = Mock()
        gene1.seq_name = 'seq1'
        self.test_gff1.genes = [gene1]
        self.assertTrue(self.test_gff1.contains_gene_on_seq('seq1'))

        
##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGFF))
    return suite

if __name__ == '__main__':
    unittest.main()

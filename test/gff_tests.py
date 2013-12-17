#!/usr/bin/env python

import unittest
from mock import Mock, patch, PropertyMock
from src.gff import GFF
from src.bed import Bed
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

    def test_validate_first_line(self):
        self.assertTrue(self.test_gff1.validate_first_line(self.test_gene_line1))
        self.assertFalse(self.test_gff1.validate_first_line(self.test_line2))

    def test_parse_attributes(self):
        expected1 = {'identifier': '2', 'name': 'BDOR_007864-RA', 'parent_id': '1'}
        string1 = 'ID=2;Name=BDOR_007864-RA;Parent=1'
        self.assertEqual(expected1, self.test_gff1.parse_attributes(string1))
        expected2 = {'identifier': '1', 'name': 'BDOR_007864'}
        string2 = 'ID=1;Name=BDOR_007864'
        self.assertEqual(expected2, self.test_gff1.parse_attributes(string2))
        bad_string1 = 'not semicolon delimited...'
        sys.stderr.write("test should generate an error message:\n")
        self.assertFalse(self.test_gff1.parse_attributes(bad_string1))
        bad_string2 = 'semicolon delimited; no equals signs though'
        sys.stderr.write("test should generate an error message:\n")
        self.assertFalse(self.test_gff1.parse_attributes(bad_string2))

    def test_extract_cds_args(self):
        expected = {'identifier': '229.1', 'name': 'BDOR_007863.1-RA:cds:44', 'indices': [10247, 10625], 'phase': 1, 'parent_id': '173.1'}
        actual = self.test_gff1.extract_cds_args(self.test_cds_line1)
        self.assertEqual(expected, actual)

    def test_extract_exon_args(self):
        expected = {'identifier': '199.1', 'name': 'BDOR_007863.1-RA:exon:48', 'indices': [10247, 10625], 'score': '100.148', 'parent_id': '173.1'}
        actual = self.test_gff1.extract_exon_args(self.test_exon_line1)
        self.assertEqual(expected, actual)

    def test_extract_other_feature_args(self):
        expected = {'feature_type': 'five_prime_UTR', 'identifier': '234.1', 'name': 'BDOR_007863.1-RA:UTR1', 'indices': [460, 627], 'parent_id': '173.1'}
        actual = self.test_gff1.extract_other_feature_args(self.test_feature_line)
        self.assertEqual(expected, actual)

    def test_extract_mrna_args(self):
        expected = {'identifier': '173.1', 'name': 'BDOR_007863.1-RA', 'indices': [460, 12713], 'parent_id': '172.1'}
        actual = self.test_gff1.extract_mrna_args(self.test_mrna_line1)
        self.assertEqual(expected, actual)

    def test_extract_gene_args(self):
        expected = {'seq_name': 'sctg_0080_0020', 'source': 'maker', 'indices': [3734, 7436], 'strand': '+', 'identifier': '1', 'name': 'BDOR_007864'}
        actual = self.test_gff1.extract_gene_args(self.test_gene_line1)
        self.assertEqual(expected, actual)

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
        # should add object to current_mrna.other_features
        self.test_gff2.process_other_feature_line(self.test_feature_line)
        self.assertEqual(1, len(self.test_gff2.current_mrna.other_features))
        self.assertEqual('five_prime_UTR', self.test_gff2.current_mrna.other_features[0].feature_type)
        

    def test_gff_file_stuff(self):
        gff = GFF()
        with open('sample_files/tiny_test.gff', 'rb') as f:
            reader = csv.reader(f, delimiter='\t')
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

    def test_apply_bed(self):
        # build mock genes
        gene1 = Mock()
        seq1 = PropertyMock(return_value = 'sctg_1')
        type(gene1).seq_name= seq1
        gene2 = Mock()
        seq2 = PropertyMock(return_value = 'sctg_2')
        type(gene2).seq_name= seq2
        # add them to the gff
        self.test_gff1.genes.append(gene1)
        self.test_gff1.genes.append(gene2)
        # easier to use a real Bed...
        bed = Bed({'sctg_1': [100, 500], 'sctg_3': [20, 50]})
        # verify that bed works as expected
        # yes i realize this isn't test_bed, but
        # could help one day with debugging?
        self.assertTrue(bed.contains('sctg_1'))
        self.assertFalse(bed.contains('sctg_2'))
        self.assertEquals([100, 500], bed.get_coordinates('sctg_1'))

        self.test_gff1.apply_bed(bed)
        assert not gene2.trim.called
        gene1.trim.assert_called_with([100, 500])

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

    def test_remove_all_gene_segments(self):
        gene1 = Mock()
        gene2 = Mock()
        gene3 = Mock()
        gene1.name = "BDOR_00001"
        gene2.name = "BDOR_00002.1"
        gene3.name = "BDOR_00002.2"
        self.assertEquals(0, len(self.test_gff1.genes))
        self.test_gff1.genes.extend([gene1, gene2, gene3])
        self.assertEquals(3, len(self.test_gff1.genes))
        self.assertEquals("BDOR_00001", self.test_gff1.genes[0].name)
        self.test_gff1.remove_all_gene_segments("BDOR_00002")
        self.assertEquals(1, len(self.test_gff1.genes))
        self.assertEquals("BDOR_00001", self.test_gff1.genes[0].name)


        
         

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGFF))
    return suite

if __name__ == '__main__':
    unittest.main()

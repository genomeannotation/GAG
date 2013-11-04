#!/usr/bin/env python

import unittest
from mock import Mock, patch
from src.gff import GFF
import sys
import csv

class TestGFF(unittest.TestCase):

    def test_gff(self):
        # test constructor
        test_gff1 = GFF()
        self.assertEqual('GFF', test_gff1.__class__.__name__)
        self.assertEqual(0, len(test_gff1.genes))
        self.assertFalse(test_gff1.current_gene)
        self.assertFalse(test_gff1.current_mrna)
        self.assertFalse(test_gff1.current_exon)
        self.assertFalse(test_gff1.current_cds)

        # test validate_line
        bad_line1 = ['sctg_x', 'maker']
        self.assertFalse(test_gff1.validate_line(bad_line1))
        test_line1 = ['sctg_0080_0020', 'maker', 'gene', '3734', '7436', '.', '+', '.', 'ID=1;Name=BDOR_007864']
        self.assertTrue(test_gff1.validate_line(test_line1))

        # test line_type
        self.assertEqual('gene', test_gff1.line_type(test_line1))
        test_line2 = ['sctg_0080_0020', 'maker', 'mRNA', '3734', '7436', '.', '+', '.', 'ID=2;Name=BDOR_007864-RA;Parent=1']
        self.assertEqual('mRNA', test_gff1.line_type(test_line2))

        # test validate_first_line
        self.assertTrue(test_gff1.validate_first_line(test_line1))
        self.assertFalse(test_gff1.validate_first_line(test_line2))

        # test parse_attributes
        expected = {'identifier': '2', 'name': 'BDOR_007864-RA', 'parent_id': '1'}
        string1 = 'ID=2;Name=BDOR_007864-RA;Parent=1'
        self.assertEqual(expected, test_gff1.parse_attributes(string1))
        expected = {'identifier': '1', 'name': 'BDOR_007864'}
        string2 = 'ID=1;Name=BDOR_007864'
        self.assertEqual(expected, test_gff1.parse_attributes(string2))
        bad_string1 = 'not semicolon delimited...'
        sys.stderr.write("test should generate an error message:\n")
        self.assertFalse(test_gff1.parse_attributes(bad_string1))
        bad_string2 = 'semicolon delimited; no equals signs though'
        sys.stderr.write("test should generate an error message:\n")
        self.assertFalse(test_gff1.parse_attributes(bad_string2))
       
        # test extract_cds_args
        test_cds_line1 = ['sctg_0080_0020', 'maker', 'CDS', '10247', '10625', '.', '-', '1', 'ID=229.1;Name=BDOR_007863.1-RA:cds:44;Parent=173.1']
        test_cds_line2 = ['sctg_0080_0020', 'maker', 'CDS', '9089', '10086', '.', '-', '2', 'ID=230.1;Name=BDOR_007863.1-RA:cds:45;Parent=173.1']
        expected = {'identifier': '229.1', 'name': 'BDOR_007863.1-RA:cds:44', 'indices': [10247, 10625], 'phase': 1, 'parent_id': '173.1'}
        actual = test_gff1.extract_cds_args(test_cds_line1)
        self.assertEqual(expected, actual)

        # test extract_exon_args
        test_exon_line1 = ['sctg_0080_0020', 'maker', 'exon', '10247', '10625', '100.148', '-', '.', 'ID=199.1;Name=BDOR_007863.1-RA:exon:48;Parent=173.1']
        test_exon_line2 = ['sctg_0080_0020', 'maker', 'exon', '9089', '10086', '270.928', '-', '.', 'ID=200.1;Name=BDOR_007863.1-RA:exon:49;Parent=173.1']
        expected = {'identifier': '199.1', 'name': 'BDOR_007863.1-RA:exon:48', 'indices': [10247, 10625], 'score': '100.148', 'parent_id': '173.1'}
        actual = test_gff1.extract_exon_args(test_exon_line1)
        self.assertEqual(expected, actual)

        # test extract_other_feature_args
        test_feature_line = ['sctg_0080_0020', 'maker', 'five_prime_UTR', '460', '627', '.', '-', '.', 'ID=234.1;Name=BDOR_007863.1-RA:UTR1;Parent=173.1']
        expected = {'feature_type': 'five_prime_UTR', 'identifier': '234.1', 'name': 'BDOR_007863.1-RA:UTR1', 'indices': [460, 627], 'parent_id': '173.1'}
        actual = test_gff1.extract_other_feature_args(test_feature_line)
        self.assertEqual(expected, actual)

        # test extract_mrna_args
        test_mrna_line1 = ['sctg_0080_0020', 'maker', 'mRNA', '460', '12713', '.', '-', '.', 'ID=173.1;Name=BDOR_007863.1-RA;Parent=172.1']
        expected = {'identifier': '173.1', 'name': 'BDOR_007863.1-RA', 'indices': [460, 12713], 'parent_id': '172.1'}
        actual = test_gff1.extract_mrna_args(test_mrna_line1)
        self.assertEqual(expected, actual)

        # test extract_gene_args
        expected = {'seq_name': 'sctg_0080_0020', 'source': 'maker', 'indices': [3734, 7436], 'strand': '+', 'identifier': '1', 'name': 'BDOR_007864'}
        actual = test_gff1.extract_gene_args(test_line1)
        self.assertEqual(expected, actual)

        # test process_line
        with patch.object(test_gff1, 'process_gene_line') as mock:
            test_gff1.process_line(test_line1)
        mock.assert_called_with(test_line1)
        with patch.object(test_gff1, 'process_mrna_line') as mock:
            test_gff1.process_line(test_mrna_line1)
        mock.assert_called_with(test_mrna_line1)
        with patch.object(test_gff1, 'process_cds_line') as mock:
            test_gff1.process_line(test_cds_line1)
        mock.assert_called_with(test_cds_line1)
        with patch.object(test_gff1, 'process_exon_line') as mock:
            test_gff1.process_line(test_exon_line1)
        mock.assert_called_with(test_exon_line1)
        with patch.object(test_gff1, 'process_other_feature_line') as mock:
            test_gff1.process_line(test_feature_line)
        mock.assert_called_with(test_feature_line)

        # test process_gene_line
        test_gff2 = GFF()
        self.assertEqual(0, len(test_gff2.genes))
        self.assertFalse(test_gff2.current_gene)
        # if no current gene, simply create new gene
        test_gff2.process_gene_line(test_line1)   
        self.assertTrue(test_gff2.current_gene)
        # if gene already present, should wrap up current
        # gene, add it to genes[] and instantiate a new one
        self.assertEqual(0, len(test_gff2.genes))
        test_gene_line2 = ['sctg_0080_0020', 'maker', 'gene', '460', '12713', '.', '-', '.', 'ID=172.1;Name=BDOR_007863.1']
        test_gff2.process_gene_line(test_gene_line2)
        self.assertEqual(1, len(test_gff2.genes))
        self.assertEqual('172.1', test_gff2.current_gene.identifier)

        # test process_mrna_line
        self.assertFalse(test_gff2.current_mrna)
        # if no current mrna, should create a new one
        test_gff2.process_mrna_line(test_mrna_line1) 
        self.assertTrue(test_gff2.current_mrna)
        # if mrna already present, should wrap up current
        # mrna, add it to current_gene, instantiate new mrna
        self.assertEqual(0, len(test_gff2.current_gene.mrnas))
        self.assertEqual('173.1', test_gff2.current_mrna.identifier)
        test_mrna_line2 = ['sctg_0080_0020', 'maker', 'mRNA', '460', '12713', '.', '-', '.', 'ID=foo;Name=BDOR_007863.1-RA;Parent=172.1']
        test_gff2.process_mrna_line(test_mrna_line2)
        self.assertEqual(1, len(test_gff2.current_gene.mrnas))
        self.assertEqual('foo', test_gff2.current_mrna.identifier)

        # test process_cds_line
        self.assertFalse(test_gff2.current_cds)
        # if no current cds should instantiate one
        test_gff2.process_cds_line(test_cds_line1)
        self.assertTrue(test_gff2.current_cds)
        self.assertEqual(['229.1'], test_gff2.current_cds.identifier)
        # if cds already present, should add info to current_cds
        test_gff2.process_cds_line(test_cds_line2)
        self.assertEqual(['229.1', '230.1'], test_gff2.current_cds.identifier)

        # test process_exon_line
        self.assertFalse(test_gff2.current_exon)
        # if no current exon should instantiate one
        test_gff2.process_exon_line(test_exon_line1)
        self.assertTrue(test_gff2.current_exon)
        self.assertEqual(1, len(test_gff2.current_exon.indices))
        # if exon already present, add info to current_exon
        test_gff2.process_exon_line(test_exon_line2)
        self.assertEqual(2, len(test_gff2.current_exon.indices))

        # test process_other_feature_line
        self.assertEqual(0, len(test_gff2.current_mrna.other_features))
        # should add object to current_mrna.other_features
        test_gff2.process_other_feature_line(test_feature_line)
        self.assertEqual(1, len(test_gff2.current_mrna.other_features))
        self.assertEqual('five_prime_UTR', test_gff2.current_mrna.other_features[0].feature_type)
        
        # REDUNDANT? but why delete a test?
        # test process_cds_line
        # should instantiate cds if no current_cds
        self.assertFalse(test_gff1.current_cds)
        test_gff1.process_cds_line(test_cds_line1)
        self.assertTrue(test_gff1.current_cds)
        # should add info to current cds if already exists
        test_gff1.current_cds = Mock()
        test_gff1.process_cds_line(test_cds_line2)
        test_gff1.current_cds.add_indices.assert_called_with([9089, 10086])
        test_gff1.current_cds.add_phase.assert_called_with(2)
        test_gff1.current_cds.add_name.assert_called_with('BDOR_007863.1-RA:cds:45')
        test_gff1.current_cds.add_identifier.assert_called_with('230.1')

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

        with open('foo', 'wb') as outfile:
            for gene in gff.genes:  
                outfile.write(gene.to_gff())
         

        # TODO verify that the above works :)  
        

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGFF))
    return suite

if __name__ == '__main__':
    unittest.main()

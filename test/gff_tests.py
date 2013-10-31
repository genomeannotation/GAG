#!/usr/bin/env python

import unittest
from mock import Mock, patch
from src.gff import GFF
import sys

class TestGFF(unittest.TestCase):

    def test_gff(self):
        # test constructor
        test_gff1 = GFF()
        self.assertEquals('GFF', test_gff1.__class__.__name__)
        self.assertEquals(0, len(test_gff1.genes))
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
        self.assertEquals('gene', test_gff1.line_type(test_line1))
        test_line2 = ['sctg_0080_0020', 'maker', 'mRNA', '3734', '7436', '.', '+', '.', 'ID=2;Name=BDOR_007864-RA;Parent=1']
        self.assertEquals('mRNA', test_gff1.line_type(test_line2))

        # test validate_first_line
        self.assertTrue(test_gff1.validate_first_line(test_line1))
        self.assertFalse(test_gff1.validate_first_line(test_line2))

        # test parse_attributes
        expected = {'identifier': '2', 'name': 'BDOR_007864-RA', 'parent_id': '1'}
        string1 = 'ID=2;Name=BDOR_007864-RA;Parent=1'
        self.assertEquals(expected, test_gff1.parse_attributes(string1))
        expected = {'identifier': '1', 'name': 'BDOR_007864'}
        string2 = 'ID=1;Name=BDOR_007864'
        self.assertEquals(expected, test_gff1.parse_attributes(string2))
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
        self.assertEquals(expected, actual)

        # test extract_exon_args
        test_exon_line1 = ['sctg_0080_0020', 'maker', 'exon', '10247', '10625', '100.148', '-', '.', 'ID=199.1;Name=BDOR_007863.1-RA:exon:48;Parent=173.1']
        test_exon_line2 = ['sctg_0080_0020', 'maker', 'exon', '9089', '10086', '270.928', '-', '.', 'ID=200.1;Name=BDOR_007863.1-RA:exon:49;Parent=173.1']
        expected = {'identifier': '199.1', 'name': 'BDOR_007863.1-RA:exon:48', 'indices': [10247, 10625], 'score': '100.148', 'parent_id': '173.1'}
        actual = test_gff1.extract_exon_args(test_exon_line1)
        self.assertEquals(expected, actual)

        # test extract_other_feature_args
        test_feature_line = ['sctg_0080_0020', 'maker', 'five_prime_UTR', '460', '627', '.', '-', '.', 'ID=234.1;Name=BDOR_007863.1-RA:UTR1;Parent=173.1']
        expected = {'feature_type': 'five_prime_UTR', 'identifier': '234.1', 'name': 'BDOR_007863.1-RA:UTR1', 'indices': [460, 627], 'parent_id': '173.1'}
        actual = test_gff1.extract_other_feature_args(test_feature_line)
        self.assertEquals(expected, actual)

        # test extract_mrna_args
        test_mrna_line1 = ['sctg_0080_0020', 'maker', 'mRNA', '460', '12713', '.', '-', '.', 'ID=173.1;Name=BDOR_007863.1-RA;Parent=172.1']
        expected = {'identifier': '173.1', 'name': 'BDOR_007863.1-RA', 'indices': [460, 12713], 'parent_id': '172.1'}
        actual = test_gff1.extract_mrna_args(test_mrna_line1)
        self.assertEquals(expected, actual)

        # test extract_gene_args
        expected = {'seq_name': 'sctg_0080_0020', 'source': 'maker', 'indices': [3734, 7436], 'strand': '+', 'identifier': '1', 'name': 'BDOR_007864'}
        actual = test_gff1.extract_gene_args(test_line1)
        self.assertEquals(expected, actual)

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
    
        
        
# read gff pseudocode:
# for each line, if comment skip
# validate line
# process line (based on line_type)
# if gene 
#   if current gene exists
#       wrap up current gene**
#   create new gene         
# if mrna
#   if current mrna exists
#       wrap up current mrna***
#   create new mrna
# if cds
#   if current cds exists
#       current_cds.update()
#   else
#       create new cds
# if exon
#   if current exon exists
#       current_exon.update()
#   else
#       create new exon
# if other feature
#   create other_feature
#   add to current mrna
# if last line of file
#   wrap up current gene

# to sum up, other_features are added to current mrna as they are created. cds and exon are added to their mrna only when new mrna or new gene is reached.
# also check parent id when add feature to mrna, or mrna to gene.

# **means: if current mrna, wrap up current mrna. add current gene to gff.genes[]. None out everything.

# ***means: if current cds/exon/etc, add to current mrna. add current mrna to current gene. None out all but gene.

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGFF))
    return suite

if __name__ == '__main__':
    unittest.main()

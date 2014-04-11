#!/usr/bin/env python

import unittest
from mock import Mock
from src.filters import *

class TestFilters(unittest.TestCase):
        
    def test_cds_length_range_filter(self):
        cds_length_range = CDSLengthRangeFilter(30, 60)
        cds_length_range.remove = False
    
        # Create a mock sequence
        seq = Mock()
        
        # Give the mock sequence some mock genes
        seq.genes = [Mock(), Mock(), Mock()]
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [Mock()]
        seq.genes[1].mrnas = [Mock(), Mock()]
        seq.genes[2].mrnas = [Mock()]
        
        # Give the mock mrnas some cds's
        seq.genes[0].mrnas[0].identifier = 'foo1-RA'
        seq.genes[0].mrnas[0].death_flagged = False
        seq.genes[0].mrnas[0].cds = Mock()
        seq.genes[0].mrnas[0].cds.length = Mock(return_value=50)
        
        seq.genes[1].mrnas[0].identifier = 'foo2-RA'
        seq.genes[1].mrnas[0].death_flagged = False
        seq.genes[1].mrnas[0].cds = None
        
        seq.genes[1].mrnas[1].identifier = 'foo2-RB'
        seq.genes[1].mrnas[1].death_flagged = False
        seq.genes[1].mrnas[1].cds = Mock()
        seq.genes[1].mrnas[1].cds.length = Mock(return_value=20)
        
        seq.genes[2].mrnas[0].identifier = 'foo3-RA'
        seq.genes[2].mrnas[0].death_flagged = False
        seq.genes[2].mrnas[0].cds = Mock()
        seq.genes[2].mrnas[0].cds.length = Mock(return_value=70)
        
        # Apply the filter
        cds_length_range.apply(seq)
        
        self.assertFalse(seq.genes[0].mrnas[0].death_flagged)
        self.assertTrue(seq.genes[1].mrnas[0].death_flagged)
        self.assertTrue(seq.genes[1].mrnas[1].death_flagged)
        self.assertTrue(seq.genes[2].mrnas[0].death_flagged)
        
        
    def test_exon_length_range_filter(self):
        exon_length_range = ExonLengthRangeFilter(30, 60)
        exon_length_range.remove = False
    
        # Create a mock sequence
        seq = Mock()
        
        # Give the mock sequence some mock genes
        seq.genes = [Mock(), Mock(), Mock()]
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [Mock()]
        seq.genes[1].mrnas = [Mock(), Mock()]
        seq.genes[2].mrnas = [Mock()]
        
        # Give the mock mrnas some exon's
        seq.genes[0].mrnas[0].identifier = 'foo1-RA'
        seq.genes[0].mrnas[0].death_flagged = False
        seq.genes[0].mrnas[0].exon = Mock()
        seq.genes[0].mrnas[0].get_shortest_exon = Mock(return_value=50)
        seq.genes[0].mrnas[0].get_longest_exon = Mock(return_value=50)
        
        seq.genes[1].mrnas[0].identifier = 'foo2-RA'
        seq.genes[1].mrnas[0].death_flagged = False
        seq.genes[1].mrnas[0].exon = Mock()
        seq.genes[1].mrnas[0].get_shortest_exon = Mock(return_value=30)
        seq.genes[1].mrnas[0].get_longest_exon = Mock(return_value=30)
        
        seq.genes[1].mrnas[1].identifier = 'foo2-RB'
        seq.genes[1].mrnas[1].death_flagged = False
        seq.genes[1].mrnas[1].exon = None
        
        seq.genes[2].mrnas[0].identifier = 'foo3-RA'
        seq.genes[2].mrnas[0].death_flagged = False
        seq.genes[2].mrnas[0].exon = Mock()
        seq.genes[2].mrnas[0].get_shortest_exon = Mock(return_value=70)
        seq.genes[2].mrnas[0].get_longest_exon = Mock(return_value=70)
        
        # Apply the filter
        exon_length_range.apply(seq)
        
        self.assertFalse(seq.genes[0].mrnas[0].death_flagged)
        self.assertFalse(seq.genes[1].mrnas[0].death_flagged)
        self.assertTrue(seq.genes[1].mrnas[1].death_flagged)
        self.assertTrue(seq.genes[2].mrnas[0].death_flagged)
        
        
    def test_intron_length_range_filter(self):
        intron_length_range = IntronLengthRangeFilter(30, 60)
        intron_length_range.remove = False
    
        # Create a mock sequence
        seq = Mock()
        
        # Give the mock sequence some mock genes
        seq.genes = [Mock(), Mock(), Mock()]
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [Mock()]
        seq.genes[1].mrnas = [Mock(), Mock()]
        seq.genes[2].mrnas = [Mock()]
        
        # Give the mock mrnas some exon's
        seq.genes[0].mrnas[0].identifier = 'foo1-RA'
        seq.genes[0].mrnas[0].death_flagged = False
        seq.genes[0].mrnas[0].exon = Mock()
        seq.genes[0].mrnas[0].get_shortest_intron = Mock(return_value=50)
        seq.genes[0].mrnas[0].get_longest_intron = Mock(return_value=50)
        
        seq.genes[1].mrnas[0].identifier = 'foo2-RA'
        seq.genes[1].mrnas[0].death_flagged = False
        seq.genes[1].mrnas[0].exon = Mock()
        seq.genes[1].mrnas[0].get_shortest_intron = Mock(return_value=30)
        seq.genes[1].mrnas[0].get_longest_intron = Mock(return_value=30)
        
        seq.genes[1].mrnas[1].identifier = 'foo2-RB'
        seq.genes[1].mrnas[1].death_flagged = False
        seq.genes[1].mrnas[1].exon = None
        
        seq.genes[2].mrnas[0].identifier = 'foo3-RA'
        seq.genes[2].mrnas[0].death_flagged = False
        seq.genes[2].mrnas[0].exon = Mock()
        seq.genes[2].mrnas[0].get_shortest_intron = Mock(return_value=70)
        seq.genes[2].mrnas[0].get_longest_intron = Mock(return_value=70)
        
        # Apply the filter
        intron_length_range.apply(seq)
        
        self.assertFalse(seq.genes[0].mrnas[0].death_flagged)
        self.assertFalse(seq.genes[1].mrnas[0].death_flagged)
        self.assertTrue(seq.genes[1].mrnas[1].death_flagged)
        self.assertTrue(seq.genes[2].mrnas[0].death_flagged)
        
    def test_gene_length_range_filter(self):
        gene_length_range = GeneLengthRangeFilter(30, 60)
        gene_length_range.remove = False
    
        # Create a mock sequence
        seq = Mock()
        
        # Give the mock sequence some mock genes
        seq.genes = [Mock(), Mock(), Mock(), Mock(), Mock()]
        
        seq.genes[0].death_flagged = False
        seq.genes[1].death_flagged = False
        seq.genes[2].death_flagged = False
        seq.genes[3].death_flagged = False
        seq.genes[4].death_flagged = False
        
        seq.genes[0].length = Mock(return_value=20)
        seq.genes[1].length = Mock(return_value=30)
        seq.genes[2].length = Mock(return_value=45)
        seq.genes[3].length = Mock(return_value=60)
        seq.genes[4].length = Mock(return_value=70)
        
        # Apply the filter
        gene_length_range.apply(seq)
        
        self.assertTrue(seq.genes[0].death_flagged)
        self.assertFalse(seq.genes[1].death_flagged)
        self.assertFalse(seq.genes[2].death_flagged)
        self.assertFalse(seq.genes[3].death_flagged)
        self.assertTrue(seq.genes[4].death_flagged)




##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFilters))
    return suite

if __name__ == '__main__':
    unittest.main()

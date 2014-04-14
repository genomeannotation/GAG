#!/usr/bin/env python

import unittest
from mock import Mock
from src.filters import *

class TestFilters(unittest.TestCase):
        
    def test_cds_length_range_filter(self):
        cds_length_range = CDSLengthRangeFilter(30, 60)
    
        # Create a mock sequence
        seq = Mock()
        
        # Give the mock sequence some mock genes
        seq.genes = [Mock(), Mock(), Mock()]
        
        # Give the mock mrnas some cds's
        test_mrna0 = Mock()
        test_mrna0.identifier = 'foo1-RA'
        test_mrna0.death_flagged = False
        test_mrna0.cds = Mock()
        test_mrna0.cds.length = Mock(return_value=50)
        
        test_mrna1 = Mock()
        test_mrna1.identifier = 'foo2-RA'
        test_mrna1.death_flagged = False
        test_mrna1.cds = None
        
        test_mrna2 = Mock()
        test_mrna2.identifier = 'foo2-RB'
        test_mrna2.death_flagged = False
        test_mrna2.cds = Mock()
        test_mrna2.cds.length = Mock(return_value=20)
        
        test_mrna3 = Mock()
        test_mrna3.identifier = 'foo3-RA'
        test_mrna3.death_flagged = False
        test_mrna3.cds = Mock()
        test_mrna3.cds.length = Mock(return_value=70)
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [test_mrna0]
        seq.genes[1].mrnas = [test_mrna1, test_mrna2]
        seq.genes[2].mrnas = [test_mrna3]
        
        # Apply the filter
        cds_length_range.apply(seq)
        
        self.assertEqual(seq.genes[0].mrnas+seq.genes[1].mrnas+seq.genes[2].mrnas, [test_mrna0])
        
        
    def test_exon_length_range_filter(self):
        exon_length_range = ExonLengthRangeFilter(30, 60)
    
        # Create a mock sequence
        seq = Mock()
        
        # Give the mock sequence some mock genes
        seq.genes = [Mock(), Mock(), Mock()]
        
        # Give the mock mrnas some exon's
        test_mrna0 = Mock()
        test_mrna0.identifier = 'foo1-RA'
        test_mrna0.death_flagged = False
        test_mrna0.exon = Mock()
        test_mrna0.get_shortest_exon = Mock(return_value=50)
        test_mrna0.get_longest_exon = Mock(return_value=50)
        
        test_mrna1 = Mock()
        test_mrna1.identifier = 'foo2-RA'
        test_mrna1.death_flagged = False
        test_mrna1.exon = Mock()
        test_mrna1.get_shortest_exon = Mock(return_value=30)
        test_mrna1.get_longest_exon = Mock(return_value=30)
        
        test_mrna2 = Mock()
        test_mrna2.identifier = 'foo2-RB'
        test_mrna2.death_flagged = False
        test_mrna2.exon = None
        
        test_mrna3 = Mock()
        test_mrna3.identifier = 'foo3-RA'
        test_mrna3.death_flagged = False
        test_mrna3.exon = Mock()
        test_mrna3.get_shortest_exon = Mock(return_value=70)
        test_mrna3.get_longest_exon = Mock(return_value=70)
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [test_mrna0]
        seq.genes[1].mrnas = [test_mrna1, test_mrna2]
        seq.genes[2].mrnas = [test_mrna3]
        
        # Apply the filter
        exon_length_range.apply(seq)
        
        self.assertEqual(seq.genes[0].mrnas+seq.genes[1].mrnas+seq.genes[2].mrnas, [test_mrna0, test_mrna1])
        
        
    def test_intron_length_range_filter(self):
        intron_length_range = IntronLengthRangeFilter(30, 60)
    
        # Create a mock sequence
        seq = Mock()
        
        # Give the mock sequence some mock genes
        seq.genes = [Mock(), Mock(), Mock()]
        
        # Give the mock mrnas some exon's
        test_mrna0 = Mock()
        test_mrna0.identifier = 'foo1-RA'
        test_mrna0.death_flagged = False
        test_mrna0.exon = Mock()
        test_mrna0.get_shortest_intron = Mock(return_value=50)
        test_mrna0.get_longest_intron = Mock(return_value=50)
        
        test_mrna1 = Mock()
        test_mrna1.identifier = 'foo2-RA'
        test_mrna1.death_flagged = False
        test_mrna1.exon = Mock()
        test_mrna1.get_shortest_intron = Mock(return_value=30)
        test_mrna1.get_longest_intron = Mock(return_value=30)
        
        test_mrna2 = Mock()
        test_mrna2.identifier = 'foo2-RB'
        test_mrna2.death_flagged = False
        test_mrna2.exon = None
        
        test_mrna3 = Mock()
        test_mrna3.identifier = 'foo3-RA'
        test_mrna3.death_flagged = False
        test_mrna3.exon = Mock()
        test_mrna3.get_shortest_intron = Mock(return_value=70)
        test_mrna3.get_longest_intron = Mock(return_value=70)
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [test_mrna0]
        seq.genes[1].mrnas = [test_mrna1, test_mrna2]
        seq.genes[2].mrnas = [test_mrna3]
        
        # Apply the filter
        intron_length_range.apply(seq)
        
        self.assertEqual(seq.genes[0].mrnas+seq.genes[1].mrnas+seq.genes[2].mrnas, [test_mrna0, test_mrna1])
        
    def test_gene_length_range_filter(self):
        gene_length_range = GeneLengthRangeFilter(30, 60)
    
        # Create a mock sequence
        seq = Mock()
        
        test_gene0 = Mock()
        test_gene1 = Mock()
        test_gene2 = Mock()
        test_gene3 = Mock()
        test_gene4 = Mock()
        
        test_gene0.death_flagged = False
        test_gene1.death_flagged = False
        test_gene2.death_flagged = False
        test_gene3.death_flagged = False
        test_gene4.death_flagged = False
        
        test_gene0.length = Mock(return_value=20)
        test_gene1.length = Mock(return_value=30)
        test_gene2.length = Mock(return_value=45)
        test_gene3.length = Mock(return_value=60)
        test_gene4.length = Mock(return_value=70)
        
        # Give the mock sequence some mock genes
        seq.genes = [test_gene0, test_gene1, test_gene2, test_gene3, test_gene4]
        
        # Apply the filter
        gene_length_range.apply(seq)
        
        self.assertEqual(seq.genes, [test_gene1, test_gene2, test_gene3])




##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFilters))
    return suite

if __name__ == '__main__':
    unittest.main()

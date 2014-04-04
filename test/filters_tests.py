#!/usr/bin/env python

import unittest
from mock import Mock
from src.cds_length_range_filter import CDSLengthRangeFilter

class TestFilters(unittest.TestCase):
        
    def test_cds_length_range_filter(self):
        cds_length_range = CDSLengthRangeFilter(30, 60)
    
        # Create a mock sequence
        seq = Mock()
        
        # Give the mock sequence some mock genes
        seq.genes = [Mock(), Mock(), Mock()]
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [Mock()]
        seq.genes[1].mrnas = [Mock(), Mock()]
        seq.genes[2].mrnas = [Mock()]
        
        # Give the mock mrnas some cds's
        seq.genes[0].mrnas[0].death_flagged = False
        seq.genes[0].mrnas[0].cds = Mock()
        seq.genes[0].mrnas[0].cds.length = Mock(return_value=50)
        
        seq.genes[1].mrnas[0].death_flagged = False
        seq.genes[1].mrnas[0].cds = None
        
        seq.genes[1].mrnas[1].death_flagged = False
        seq.genes[1].mrnas[1].cds = Mock()
        seq.genes[1].mrnas[1].cds.length = Mock(return_value=20)
        
        seq.genes[2].mrnas[0].death_flagged = False
        seq.genes[2].mrnas[0].cds = Mock()
        seq.genes[2].mrnas[0].cds.length = Mock(return_value=70)
        
        # Apply the filter
        cds_length_range.apply(seq)
        
        self.assertFalse(seq.genes[0].mrnas[0].death_flagged)
        self.assertFalse(seq.genes[1].mrnas[0].death_flagged)
        self.assertTrue(seq.genes[1].mrnas[1].death_flagged)
        self.assertTrue(seq.genes[2].mrnas[0].death_flagged)




##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFilters))
    return suite

if __name__ == '__main__':
    unittest.main()

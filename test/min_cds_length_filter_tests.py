#!/usr/bin/env python

import unittest
from mock import Mock
from src.min_cds_length_filter import MinCDSLengthFilter

class TestSequence(unittest.TestCase):

    def setUp(self):
        self.filter = MinCDSLengthFilter(30)
        
    def test_filter(self):
        # Create a mock sequence
        seq = Mock()
        
        # Give the mock sequence some mock genes
        seq.genes = [Mock(), Mock()]
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [Mock()]
        seq.genes[1].mrnas = [Mock(), Mock()]
        
        # Give the mock mrnas some cds's
        seq.genes[0].mrnas[0].cds = Mock()
        seq.genes[0].mrnas[0].cds.length = Mock(return_value=50)
        
        seq.genes[1].mrnas[0].cds = None
        
        seq.genes[1].mrnas[1].cds = Mock()
        seq.genes[1].mrnas[1].cds.length = Mock(return_value=20)
        
        # Apply the filter
        self.filter.apply(seq)
        
        self.assertNotEqual(seq.genes[0].mrnas[0].cds, None)
        self.assertEqual(seq.genes[1].mrnas[1].cds, None)




##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSequence))
    return suite

if __name__ == '__main__':
    unittest.main()

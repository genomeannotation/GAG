#!/usr/bin/env python

import unittest
from mock import Mock
from src.bed import Bed
import sys
import csv

class TestBed(unittest.TestCase):

    def test_bed(self):
        # test constructor
        test_bed1 = Bed()
        self.assertEqual('Bed', test_bed1.__class__.__name__)
        self.assertFalse(test_bed1.entries)

        # test .add_entry
        test_key1 = "sctg_0002_0347"
        test_vals1 = [1, 4307]
        test_bed1.add_entry(test_key1, test_vals1)
        self.assertEqual(1, test_bed1.entries['sctg_0002_0347'][0])
        # what if bad input?
        self.assertRaises(TypeError, test_bed1.add_entry, ['foo', 'bar'])
        self.assertRaises(TypeError, test_bed1.add_entry, ['foo', [1, 'f']])

        # test .contains
        self.assertTrue(test_bed1.contains('sctg_0002_0347'))
        self.assertFalse(test_bed1.contains('totally made-up key'))

        # test .get_coordinates
        self.assertEqual([1, 4307], test_bed1.get_coordinates('sctg_0002_0347'))
        # what if bad input?
        self.assertEqual(None, test_bed1.get_coordinates('no such seq'))

        # test .process_line
        self.assertEqual(1, len(test_bed1.entries))
        test_input1 = ['sctg_0002_0348', '1', '102']
        test_bed1.process_line(test_input1)
        self.assertEqual(2, len(test_bed1.entries))
        # what if bad input?
        bad_input1 = ['sctg_0002_0348', '7']
        sys.stderr.write("test should generate an error message:\n")
        test_bed1.process_line(bad_input1)
        self.assertEqual(2, len(test_bed1.entries))

        # test .read_file
        test_bed2 = Bed()
        with open('sample_files/sample1.bed', 'rb') as f:
            reader = csv.reader(f, delimiter='\t')
            test_bed2.read_file(reader)
        self.assertEqual(2, len(test_bed2.entries))
        self.assertTrue(test_bed2.contains('sctg_0002_0028'))
        self.assertTrue(test_bed2.contains('sctg_0034_0004'))
        coords = test_bed2.get_coordinates('sctg_0034_0004')
        self.assertEqual(4431, coords[1])
        

        
        



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestBed))
    return suite

if __name__ == '__main__':
    unittest.main()

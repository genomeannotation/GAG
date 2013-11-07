#!/usr/bin/env python

import unittest
from mock import Mock
from src.bed import Bed
import sys
import csv

class TestBed(unittest.TestCase):
    def setUp(self):
        self.test_bed0 = Bed()
        self.test_bed1 = Bed()
        self.test_bed2 = Bed()
        test_key1 = "sctg_0002_0347"
        test_vals1 = [1, 4307]
        self.test_bed1.add_entry(test_key1, test_vals1)

    def test_constructor(self):
        self.assertEqual('Bed', self.test_bed0.__class__.__name__)
        self.assertFalse(self.test_bed0.entries)

    def test_add_entry(self):
        test_key1 = "sctg_0002_0347"
        test_vals1 = [1, 4307]
        self.test_bed0.add_entry(test_key1, test_vals1)
        self.assertEqual(1, self.test_bed0.entries['sctg_0002_0347'][0])
        # what if bad input?
        bad_input1 = ['foo', 'bar']
        bad_input2 = ['foo', [1, 'f']]
        self.assertRaises(TypeError, self.test_bed0.add_entry, bad_input1)
        self.assertRaises(TypeError, self.test_bed0.add_entry, bad_input2)

    def test_contains(self):
        self.assertTrue(self.test_bed1.contains('sctg_0002_0347'))
        self.assertFalse(self.test_bed1.contains('totally made-up key'))

    def test_get_coordinates(self):
        self.assertEqual([1, 4307], self.test_bed1.get_coordinates('sctg_0002_0347'))
        # what if bad input?
        self.assertEqual(None, self.test_bed1.get_coordinates('no such seq'))

    def test_process_line(self):
        self.assertEqual(1, len(self.test_bed1.entries))
        test_input1 = ['sctg_0002_0348', '1', '102']
        self.test_bed1.process_line(test_input1)
        self.assertEqual(2, len(self.test_bed1.entries))
        # what if bad input?
        bad_input1 = ['sctg_0002_0348', '7']
        sys.stderr.write("test should generate an error message:\n")
        self.test_bed1.process_line(bad_input1)
        self.assertEqual(2, len(self.test_bed1.entries))

    def test_read_file(self):
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

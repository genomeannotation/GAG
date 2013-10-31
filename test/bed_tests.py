#!/usr/bin/env python

import unittest
from mock import Mock
from src.bed import Bed

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



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestBed))
    return suite

if __name__ == '__main__':
    unittest.main()

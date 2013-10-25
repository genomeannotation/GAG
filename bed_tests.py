#!/usr/bin/env python

import unittest
from mock import Mock
from bed import Bed

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

        # test .contains
        self.assertTrue(test_bed1.contains('sctg_0002_0347'))
        self.assertFalse(test_bed1.contains('totally made-up key'))



##########################
if __name__ == '__main__':
    unittest.main()

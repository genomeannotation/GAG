#!/usr/bin/env python

import unittest
from gff_reader import GffReader

class TestStuff(unittest.TestCase):


    def test_gff_reader(self):
      test_reader = GffReader()
      stats = test_reader.summary_stats("test_files/test.gff")
      self.assertEqual("test_files/test.gff: 14 lines, 1 genes, 1 mRNA, 5 exons, 5 CDS, 1 start codons, 1 stop codons", stats)


##########################
if __name__ == '__main__':
    unittest.main()

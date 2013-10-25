#!/usr/bin/env python

import unittest
from gff_reader import GffReader
from fasta_reader import FastaReader
from trinotate_reader import TrinotateReader
from feature_tbl_writer import FeatureTblWriter
from feature_tbl_entry import FeatureTblEntry
from fasta import Fasta
import sqlite3

class TestStuff(unittest.TestCase):


    def test_fasta(self):
        fasta = Fasta()
        fasta.readString('>seq1\nATGCCGTA\n>seq2\nAGGTCC\n>seq3\nGGGGGG')
        fasta.trimSeq('seq2', 2, 3)
        self.assertEqual(fasta.writeString(), '>seq1\nATGCCGTA\n>seq2\nATCC\n>seq3\nGGGGGG')
        fasta.trimSeq('seq2', 1, 4)
        self.assertEqual(fasta.writeString(), '>seq1\nATGCCGTA\n>seq3\nGGGGGG')


##########################
if __name__ == '__main__':
    unittest.main()

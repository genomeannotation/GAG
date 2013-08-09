#!/usr/bin/env python

import unittest
from gff_reader import GffReader
from fasta_reader import FastaReader
from sqlite_wrapper import SqliteWrapper
import sqlite3

class TestStuff(unittest.TestCase):


    def test_sqlite_wrapper(self):
		db = sqlite3.connect(':memory:')
		sqlite = SqliteWrapper(db)
		sqlite.createTable('sample', ['id INTEGER PRIMARY KEY', 'name TEXT'])
		sqlite.insertRow('sample', [1, '"hello"'])
		row = sqlite.getRowAsStr('sample', 'id', '1')
		self.assertEqual(row, "(1, u'hello')")


    def test_gff_reader(self):
        test_reader = GffReader()
        stats = test_reader.summary_stats("test_files/test.gff")
        self.assertEqual("test_files/test.gff: 14 lines, 1 genes, 1 mRNA, 5 exons, 5 CDS, 1 start codons, 1 stop codons", stats)


    def test_gff_reader2(self):
        test_reader = GffReader()
        db = test_reader.load('test_files/test.gff')
        c = db.cursor()
        c.execute('SELECT * FROM gff')


    def test_fasta_reader(self):
        test_reader = FastaReader()
        test_reader.read_sequences_into_db("test_files/test.fasta", "test.db")
        con = sqlite3.connect("test.db")
        cur = con.cursor()
        cur.execute("select * from sequence")
        first_row = cur.fetchone()
        first_sequence_id = first_row[0]
        self.assertEqual("scaffold00080", first_sequence_id)
        first_sequence = first_row[1]
        self.assertEqual("ggaattCAgg", first_sequence[0:10])
        second_row = cur.fetchone()
        second_sequence_id = second_row[0]
        self.assertEqual("scaffold00081", second_sequence_id)
        second_sequence = second_row[1]
        self.assertEqual("gcagtttttCGCcGAAAACCCAGAAAAATGGCAAG", second_sequence[0:35])


##########################
if __name__ == '__main__':
    unittest.main()

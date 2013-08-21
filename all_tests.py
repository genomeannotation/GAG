#!/usr/bin/env python

import unittest
from gff_reader import GffReader
from fasta_reader import FastaReader
from trinotate_reader import TrinotateReader
from feature_tbl_writer import FeatureTblWriter
from sqlite_wrapper import SqliteWrapper
import sqlite3

class TestStuff(unittest.TestCase):


    def test_sqlite_wrapper(self):
        sqlite = SqliteWrapper(':memory:')
        sqlite.createTable('sample', 'id INTEGER PRIMARY KEY, name TEXT')
        sqlite.insertRow('sample', [1, "hello"])
        row = sqlite.getRowAsStr('sample', 'id', '1')
        self.assertEqual(row, "(1, u'hello')")


    def test_gff_reader(self):
        test_reader = GffReader()
        stats = test_reader.summary_stats("test_files/test.gff")
        self.assertEqual("test_files/test.gff: 14 lines, 1 genes, 1 mRNA, 5 exons, 5 CDS, 1 start codons, 1 stop codons", stats)


    def test_gff_reader2(self):
        gff_db = sqlite3.connect(':memory:')
        test_reader = GffReader()
        test_reader.read_into_db('test_files/test.gff', gff_db)
        c = gff_db.cursor()
        c.execute('SELECT * FROM gff WHERE id=5')
        row = c.fetchone()
        self.assertEqual(row, ('5', 'scaffold00080', 'maker', 'exon', 106816, 107602, '0.9', '+', '.', 'BDOR_007864-RA:exon:2', '2'))


    def test_fasta_reader(self):
        con = sqlite3.connect(':memory:')
        test_reader = FastaReader()
        test_reader.read_into_db("test_files/test.fasta", con)
        cur = con.cursor()
        cur.execute("select * from fasta")
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

    def test_trinotate_reader(self):
        con = sqlite3.connect(':memory:')
        test_reader = TrinotateReader()
        test_reader.read_into_db("test_files/test.trinotate.xls", con)
        cur = con.cursor()
        cur.execute("select * from trinotate")
        first_row = cur.fetchone()
        self.assertEqual('BDOR_007864-RA', first_row[0])
        self.assertEqual('comp10010_c0', first_row[1])
        self.assertEqual('comp10010_c0_seq1:3-1250(-)', first_row[2])
        self.assertEqual('sp|Q7Z4T9|AAT1_HUMAN`Q7Z4T9`Q:533-644,H:400-513`37.39%ID`E:9e-08`RecName: Full=Protein MAATS1; AltName: Full=AMY-1-associating protein expressed in testis 1; Short=AAT-1; AltName: Full=MYCBP/AMY-1-associated testis-expressed protein 1;`Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo.', first_row[3])
        self.assertEqual('.', first_row[4])
        self.assertEqual('.', first_row[5])
        self.assertEqual('.', first_row[6])
        self.assertEqual('NOG84773^ Amy-1-Associating protein expressed in testis 1 (AAT-1)', first_row[7])
        self.assertEqual('GO:0005739^cellular_component^mitochondrion', first_row[8])
        self.assertEqual('MAQPSHHKLAKRVFKSVYDAECSTADEKSDLKMFKNARNKATTAK', first_row[9][:45])

    def test_tbl_writer(self):
        con = sqlite3.connect(':memory:')

        gff_reader = GffReader()
        gff_reader.read_into_db('test_files/test.gff', con)

        fasta_reader = FastaReader()
        fasta_reader.read_into_db("test_files/test.fasta", con)
 
        test_writer = FeatureTblWriter()
        test_writer.write_to_db(con)
        c = con.cursor()
        c.execute('SELECT * FROM tbl WHERE type="gene"')
        row = c.fetchone()
        self.assertEqual(row, ('BDOR_007864-RA', 'scaffold00080', 'gene', '106151', '109853', '+', '.'))


##########################
if __name__ == '__main__':
    unittest.main()

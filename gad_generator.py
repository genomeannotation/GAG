#!/usr/bin/env python

import time
import sqlite3
from gff_reader import GffReader
from fasta_reader import FastaReader
from trinotate_reader import TrinotateReader
from feature_tbl_writer import FeatureTblWriter
import os

start_time = time.time()

os.system('rm real_files/tbl_db.sqlite')
con = sqlite3.connect('real_files/tbl_db.sqlite')

print("Reading gff...")
gff_reader = GffReader()
gff_reader.read_into_db('real_files/real_output.gff', con)

print(time.time() - start_time, "seconds")

print("Reading fasta...")
fasta_reader = FastaReader()
fasta_reader.read_into_db("real_files/454ScaffoldContigs.fna", con)

print(time.time() - start_time, "seconds")

print("Reading trinotate...")
trinotate_reader = TrinotateReader()
trinotate_reader.read_into_db("real_files/maker.xls", con)

print(time.time() - start_time, "seconds")

print("Writing tbl database...")
test_writer = FeatureTblWriter()
#test_writer.write_to_db(con)
test_writer.write_to_file(con, 'bdor_genome.tbl')

con.commit()

#cur = con.cursor()
#cur.execute('SELECT * FROM tbl')
#for row in cur.fetchall():
#    print(row)

print(time.time() - start_time, "seconds")

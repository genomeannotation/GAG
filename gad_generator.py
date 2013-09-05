#!/usr/bin/env python

import time
import sqlite3
from gff_reader import GffReader
from fasta_reader import FastaReader
from trinotate_reader import TrinotateReader
from feature_tbl_writer import FeatureTblWriter
import os
import sys

usage_message = """Usage: python gad_generator.py <project_name> \n\
where the working directory contains files\n\
project_name.gff, project_name.fasta and \
project_name.xls"""

if len(sys.argv) != 2:
    print(usage_message)
    sys.exit()

project_name = sys.argv[1]
sqlite_database = project_name + '.sqlite'
gff_file = project_name + '.gff'
fasta_file = project_name + '.fasta'
trinotate_file = project_name + '.xls'
tbl_file = project_name + '.tbl'

start_time = time.time()

if os.path.isfile(sqlite_database):
    os.system('rm '+sqlite_database)

con = sqlite3.connect(sqlite_database)

print("Reading gff...")
gff_reader = GffReader()
gff_reader.read_into_db(gff_file, con)

print(time.time() - start_time, "seconds")

print("Reading fasta...")
fasta_reader = FastaReader()
fasta_reader.read_into_db(fasta_file, con)

print(time.time() - start_time, "seconds")

print("Reading trinotate...")
trinotate_reader = TrinotateReader()
trinotate_reader.read_into_db(trinotate_file, con)

print(time.time() - start_time, "seconds")

print("Writing tbl database...")
test_writer = FeatureTblWriter()
#test_writer.write_to_db(con)
test_writer.write_to_file(con, tbl_file )

con.commit()

#cur = con.cursor()
#cur.execute('SELECT * FROM tbl')
#for row in cur.fetchall():
#    print(row)

print(time.time() - start_time, "seconds")

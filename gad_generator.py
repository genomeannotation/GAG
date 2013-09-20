#!/usr/bin/env python

import time
import sqlite3
from gff_reader import GffReader
from fasta_reader import FastaReader
from trinotate_reader import TrinotateReader
from feature_tbl_writer import FeatureTblWriter
import os
import sys
import glob
from sets import Set

usage_message = """Usage: python gad_generator.py <directory>\n\
where directory contains files with extensions\n\
.gff, .fasta and .xls"""

def get_file_with_extension(ext, opt = false):
    results = glob.glob('./*.' + ext)
    if len(results) != 1:
        print("Error -- directory " + working_dir + " should contain exactly one file of type " + ext)
        sys.exit()
    else:
        return results[0]

# Validate command line argument
if len(sys.argv) < 2:
    print(usage_message)
    sys.exit()

working_dir = sys.argv[1]
if not os.path.isdir(working_dir):
    print("Sorry, couldn't find directory " + working_dir)            

os.chdir(working_dir)
gff_file = get_file_with_extension('gff')
fasta_file = get_file_with_extension('fasta')
trinotate_file = get_file_with_extension('xls')
blacklist_file = get_file_with_extension('blacklist')
sqlite_database = 'GAD.sqlite'
tbl_file = 'GAD_generator_output.tbl'

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

# Get the blacklist
gene_blacklist = Set()
for badgene in blacklist_file:
    gene_blacklist.add(badgene)

print("Writing tbl database...")
test_writer = FeatureTblWriter()
#test_writer.write_to_db(con)
test_writer.write_to_file(con, tbl_file, gene_blacklist)

con.commit()

#cur = con.cursor()
#cur.execute('SELECT * FROM tbl')
#for row in cur.fetchall():
#    print(row)

print(time.time() - start_time, "seconds")

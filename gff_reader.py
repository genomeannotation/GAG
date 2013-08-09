#!/usr/bin/env python

import csv
import sqlite3

class GffReader:

	## Loads the gff file into a sqlite database
	def load(self, filename):
		database = sqlite3.connect(':memory:')
		c = database.cursor()

		# Create the sqlite database: id | seq_id | source | type | start | stop | score | strand | phase | name | parent
		c.execute('CREATE TABLE gff(id TEXT PRIMARY KEY, seq_id TEXT, source TEXT, type TEXT, start INTEGER, stop INTEGER, score TEXT, strand TEXT, phase TEXT, name TEXT, parent TEXT)')

		with open(filename, 'rb') as gff:
			reader = csv.reader(gff, delimiter='\t', quotechar='|')
			entry_id = 0
			for line in reader:
				entry_id = entry_id+1
				entry_name = 'hello'
				entry_parent = '1'
				c.execute('INSERT INTO gff VALUES("'+str(entry_id)+'","'+line[0]+'","'+line[1]+'","'+line[2]+'",'+line[3]+','+line[4]+',"'+line[5]+'","'+line[6]+'","'+line[7]+'","'+entry_name+'","'+entry_parent+'")')
		return database

    ## mostly useless function, just here to demonstrate gff-reading and unit test setup...
	def summary_stats(self, filename):
		line_count = 0
		gene_count = 0
		mRNA_count = 0
		exon_count = 0
		CDS_count = 0
		start_codon_count = 0
		stop_codon_count = 0
		with open(filename, 'rb') as gff:
			reader = csv.reader(gff, delimiter='	', quotechar='|')
			for line in reader:
				line_count += 1
				feature_type = line[2]
				if feature_type == "gene":
					gene_count += 1
				elif feature_type == "mRNA":
					mRNA_count += 1
				elif feature_type == "exon":
					exon_count += 1
				elif feature_type == "CDS":
					CDS_count += 1
				elif feature_type == "start_codon":
					start_codon_count += 1
				elif feature_type == "stop_codon":
					stop_codon_count += 1
				else:
					sys.stderr.write("Warning: unknown feature type")
		return (filename + ": " + str(line_count) + " lines, " + str(gene_count) + " genes, " + str(mRNA_count) + " mRNA, " +
					str(exon_count) + " exons, " + str(CDS_count) + " CDS, " + str(start_codon_count) + " start codons, " + 
					str(stop_codon_count) + " stop codons")

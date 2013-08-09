#!/usr/bin/env python

import sqlite3
import Bio
from Bio import SeqIO

class FastaReader:


	def read_sequences_into_db(self, filename, db_name):
		conn = sqlite3.connect(db_name)
		c = conn.cursor()
		c.execute('''CREATE TABLE IF NOT EXISTS sequence
					(seq_id text, sequence text)''')

		for seq_record in SeqIO.parse(filename, "fasta"):
			id_and_sequence = [str(seq_record.id), str(seq_record.seq)]
			c.execute('''INSERT INTO sequence VALUES (?, ?)''', id_and_sequence)
			conn.commit()

#!/usr/bin/env python

from sqlite_wrapper import SqliteWrapper

class FastaReader:


	def read_sequences_into_db(self, filename, db_name):
		# Create the sqlite database: seq_id | sequence
		sqlite = SqliteWrapper(db_name)
		sqlite.createTable('fasta', 'seq_id TEXT PRIMARY KEY, sequence TEXT')
		
		f = open(filename, 'r')
		while True:
			# Get the >
			c = f.read(1)
			if c != '>':
				raise IOError('Invalid fasta file.')

			# Get the seq_id
			seq_id = f.readline().strip().split()[0]

			# Get the sequence
			seq = str()
			last_pos = f.tell()
			c = f.read(1)
			while c != '>' and len(c) == 1:
				seq = seq+c
				last_pos = f.tell()
				c = f.read(1)

			# Write to the database
			sqlite.insertRow('fasta', [seq_id, seq])

			if len(c) < 1:
				break # Reached end of file
			else:
				f.seek(last_pos) # Put the > back so we can read the next sequence
		sqlite.commit()
		return sqlite
		
		

		

#!/usr/bin/env python

import sqlite3

class FastaReader:


    def read_sequences_into_db(self, filename, db_name):
        # Create the sqlite database: seq_id | sequence
        db_conn = sqlite3.connect(db_name)
        db_cur = db_conn.cursor()
        db_cur.execute('CREATE TABLE fasta(seq_id TEXT PRIMARY KEY, sequence TEXT)')

        with open(filename, 'r') as f:
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
                db_cur.execute('INSERT INTO fasta VALUES(?, ?)', [seq_id, seq])
                #sqlite.insertRow('fasta', [seq_id, seq])

                if len(c) < 1:
                    break # Reached end of file
                else:
                    f.seek(last_pos) # Put the > back so we can read the next sequence
            f.close()
        return db_conn
		
		

		

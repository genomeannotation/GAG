#!/usr/bin/env python

import sqlite3

class FastaReader:


    def read_into_db(self, file_name, db_conn):
        # Create the sqlite database: seq_id | sequence
        db_cur = db_conn.cursor()
        db_cur.execute('CREATE TABLE fasta(seq_id TEXT PRIMARY KEY, seq TEXT)')

        with open(file_name, 'r') as f:
            seq_id = ''
            seq = ''
            for line in f:
                if line[0] == '>':
                    if len(seq_id) > 0:
                        # Write to the database
                        db_cur.execute('INSERT INTO fasta VALUES(?, ?)', [seq_id, seq])
                    
                    seq_id = line[1:].strip().split()[0] # Get the next seq_id
                    seq = ''
                else:
                    seq += line.strip()
            # Add the last sequence to the database
            db_cur.execute('INSERT INTO fasta VALUES(?, ?)', [seq_id, seq])

        return db_conn
		
		

		

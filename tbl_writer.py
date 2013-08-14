#!/usr/bin/env python

import sqlite3

class TblWriter:

    def write_to_db(self, db_conn):

        # Get the database cursor
        db_cur = db_conn.cursor()

        # Create a new table for our tbl file
        db_cur.execute('CREATE TABLE tbl(gene_id TEXT PRIMARY KEY, seq_id TEXT, starts TEXT, stops TEXT)')

        # Fetch info from the gff table
        try:
            db_cur.execute('SELECT * FROM gff')
        except OperationalError:
            print("ERROR: gff table doesn't exist in database "+db_name)
            return

        rows = db_cur.fetchall()
        for row in rows:
            db_cur.execute('SELECT EXISTS(SELECT 1 FROM tbl WHERE gene_id=? LIMIT 1)', [row[9]])
            exists = db_cur.fetchone()[0]
            if exists == 0:
                db_cur.execute('INSERT INTO tbl VALUES(?, ?, ?, ?)', [row[9], row[1], row[4], row[5]])

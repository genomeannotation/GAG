#!/usr/bin/env python

import sqlite3

class FastaFilter:

    def __init__(self):
        self.seqStart = 0
        self.seqStop = 0

    def setIndices(self, start, stop):
        self.seqStart = start
        self.seqStop = stop

    def doQuery(self, db_conn):
        db_cur = db_conn.cursor()
        
        query = 'SELECT seq_id, substr(seq, '
        query += str(self.seqStart)+', '
        if self.seqStop >= 0:
            query += str(self.seqStop-self.seqStart)
        else:
            query += '(len(seq)-1+'+str(self.seqStop)+')-'+str(self.seqStart)
        query += ')'
                
        db_cur.execute(query)
        return db_cur.fetchall()
		
		

		

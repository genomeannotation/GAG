#!/usr/bin/env python

import csv
import sqlite3

class TrinotateReader:

    ## Loads the gff file into a sqlite database
    def read_into_db(self, file_name, db_conn):
        # Create the sqlite database: component | trans_derived | prot_id | top_blast_hit, pfam, tmhmm, eggnog, gene_ontology, prot_seq
        db_cur = db_conn.cursor()
        db_cur.execute('CREATE TABLE trinotate(prot_id TEXT PRIMARY KEY, component TEXT, trans_derived TEXT, top_blast_hit TEXT, pfam TEXT, tmhmm TEXT, eggnog TEXT, gene_ontology TEXT, prot_seq TEXT)')

        with open(file_name, 'r') as gff:
            reader = csv.reader(gff, delimiter='\t')
            for line in reader:
                db_cur.execute('INSERT INTO trinotate VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', [line[2], line[0], line[1], line[3], line[4], line[5], line[6], line[7], line[8], line[9]])      

        return db_conn

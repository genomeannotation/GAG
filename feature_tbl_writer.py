#!/usr/bin/env python

import sqlite3
import re

def regexp(expr, item):
    reg = re.compile(expr)
    return reg.search(item) is not None

class FeatureTblWriter:

    def write_to_db(self, db_conn):

        db_conn.create_function('regexp', 2, regexp)

        # Get the database cursor
        db_cur = db_conn.cursor()

        # Create a new table for our tbl file
        db_cur.execute('CREATE TABLE tbl(prot_id TEXT, seq_id TEXT, type TEXT, starts TEXT, stops TEXT, has_start INT, has_stop INT, strand TEXT, frame INT)')

        # Select all of the mRNA entries from the gff table
        try:
            db_cur.execute('SELECT * FROM gff WHERE type="mRNA"')
        except OperationalError:
            print("ERROR: gff table doesn't exist in database "+db_name)
            return

        # Iterate through the mRNA
        rnas = db_cur.fetchall()
        for rna in rnas:
            db_cur.execute('SELECT * FROM gff WHERE type="gene" AND id=? LIMIT 1', [rna[10]])
            gene = db_cur.fetchone()

            db_cur.execute('SELECT EXISTS(SELECT prot_id FROM tbl WHERE prot_id=? LIMIT 1)', [rna[9]]) # Make sure this gene entry doesn't already exist
            exists = db_cur.fetchone()[0]
            if exists == 1:
                print("ERROR: Gene already exists: ", rna[9])
                continue

            # Get whether the sequence has start and stop codons
            db_cur.execute('SELECT EXISTS(SELECT type, name FROM gff WHERE type="start_codon" AND name REGEXP ? LIMIT 1)', [rna[9]])
            has_start = db_cur.fetchone()[0]
            db_cur.execute('SELECT EXISTS(SELECT type, name FROM gff WHERE type="stop_codon" AND name REGEXP ? LIMIT 1)', [rna[9]])
            has_stop = db_cur.fetchone()[0]

            # Add the gene entry to our table
            db_cur.execute('INSERT INTO tbl VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)', [rna[9], gene[1], 'gene', gene[4], gene[5], has_start, has_stop, gene[7], gene[8]])

            # Fetch all of the exons pertaining to this gene
            db_cur.execute('SELECT * FROM gff WHERE type="exon" AND name REGEXP ?', [rna[9]])
            exons = db_cur.fetchall()
            exon_starts = ''
            exon_stops = ''
            for exon in exons:
                if len(exon_starts) > 0:
                    exon_starts += ','
                    exon_stops += ','
                exon_starts += str(exon[4])
                exon_stops += str(exon[5])

            # Add the exon entry to our table
            db_cur.execute('INSERT INTO tbl VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)', [rna[9], gene[1], 'exon', exon_starts, exon_stops, has_start, has_stop, gene[7], gene[8]])

            # Fetch all of the CDS pertaining to this gene
            db_cur.execute('SELECT * FROM gff WHERE type="CDS" AND name REGEXP ?', [rna[9]])
            cdss = db_cur.fetchall()
            cds_starts = ''
            cds_stops = ''
            for cds in cdss:
                if len(cds_starts) > 0:
                    cds_starts += ','
                    cds_stops += ','
                cds_starts += str(cds[4])
                cds_stops += str(cds[5])

            # Add the CDS entry to our table
            db_cur.execute('INSERT INTO tbl VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)', [rna[9], gene[1], 'CDS', cds_starts, cds_stops, has_start, has_stop, gene[7], gene[8]])

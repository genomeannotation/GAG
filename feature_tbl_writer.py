#!/usr/bin/env python

import sqlite3
import re

def regexp(expr, item):
    reg = re.compile(expr)
    return reg.search(item) is not None

# Takes gene ontology and returns it as a CSV of key=value
def parse_gene_ontology(gene_ontology):
    if len(gene_ontology) == 0 or gene_ontology == '.':
        return ''

    attributes = ''
    elements = gene_ontology.split('`')
    for element in elements:
        if len(attributes) > 0:
            attributes += ','

        vals = element.split('^')
        if vals[1] == 'molecular_function':
            attributes += 'go_function='
        elif vals[1] == 'cellular_component':
            attributes += 'go_component='
        elif vals[1] == 'biological_process':
            attributes += 'go_process='
        attributes += vals[2]+'|'+vals[0][3:]+'||IEA'
    return attributes

class FeatureTblWriter:

    def write_to_db(self, db_conn):

        db_conn.create_function('regexp', 2, regexp)

        # Get the database cursor
        db_cur = db_conn.cursor()

        # Create a new table for our tbl file
        db_cur.execute('CREATE TABLE tbl(prot_id TEXT, seq_id TEXT, type TEXT, starts TEXT, stops TEXT, has_start INT, has_stop INT, strand TEXT, frame INT, annotations TEXT)')

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
            
            # Get the trinotate stuff
            db_cur.execute('SELECT * FROM trinotate WHERE prot_id=? LIMIT 1', [rna[9]])
            trinotate = db_cur.fetchone()

            # Get whether the sequence has start and stop codons
            db_cur.execute('SELECT EXISTS(SELECT type, name FROM gff WHERE type="start_codon" AND name REGEXP ? LIMIT 1)', [rna[9]])
            has_start = db_cur.fetchone()[0]
            db_cur.execute('SELECT EXISTS(SELECT type, name FROM gff WHERE type="stop_codon" AND name REGEXP ? LIMIT 1)', [rna[9]])
            has_stop = db_cur.fetchone()[0]

            # Add the gene entry to our table
            db_cur.execute('INSERT INTO tbl VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', [rna[9], gene[1], 'gene', gene[4], gene[5], has_start, has_stop, gene[7], gene[8], ''])

            # Fetch all of the exons pertaining to this gene
            db_cur.execute('SELECT * FROM gff WHERE type="exon" AND name REGEXP ?', [rna[9]])
            exons = db_cur.fetchall()
            if len(exons) > 0: # Make sure there is at least one exon entry
                exon_starts = ''
                exon_stops = ''
                for exon in exons:
                    if len(exon_starts) > 0:
                        exon_starts += ','
                        exon_stops += ','
                    exon_starts += str(exon[4])
                    exon_stops += str(exon[5])

                # Add the exon entry to our table
                db_cur.execute('INSERT INTO tbl VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', [rna[9], gene[1], 'exon', exon_starts, exon_stops, has_start, has_stop, gene[7], gene[8], ''])

            # Fetch all of the CDS pertaining to this gene
            db_cur.execute('SELECT * FROM gff WHERE type="CDS" AND name REGEXP ?', [rna[9]])
            cdss = db_cur.fetchall()
            if len(cdss) > 0: # Make sure there is at least one cds entry
                cds_starts = ''
                cds_stops = ''
                for cds in cdss:
                    if len(cds_starts) > 0:
                        cds_starts += ','
                        cds_stops += ','
                    cds_starts += str(cds[4])
                    cds_stops += str(cds[5])

                # Add the CDS entry to our table
                cds_ann = parse_gene_ontology(trinotate[8])
                db_cur.execute('INSERT INTO tbl VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', [rna[9], gene[1], 'CDS', cds_starts, cds_stops, has_start, has_stop, gene[7], gene[8], cds_ann])



    def write_to_file(self, db_conn, fileName):
        db_cur = db_conn.cursor()

        with open(fileName, 'w') as f:
            # Get the sequence ids
            db_cur.execute('SELECT seq_id FROM fasta')
            seq_ids = []
            for seq_id in db_cur.fetchall():
                seq_ids.append(seq_id[0])

            print("Sequences: "+str(len(seq_ids)))
            i = 0

            for seq_id in seq_ids:
                i += 1
                if i%100 == 0:
                    print(i)

                f.write('>'+seq_id+'\n')

                # Create the temporary table
                # db_cur.execute('CREATE TABLE tmp_cur_seq(id TEXT PRIMARY KEY, seq_id TEXT, type TEXT, start INT, stop INT, name TEXT, parent TEXT)')
                
                magic_query = "CREATE TEMP TABLE tmp_cur_seq AS SELECT id, seq_id, type, start, stop, name, parent FROM gff WHERE seq_id=?"

                # Store temporary table of all info pertaining to this sequence
                db_cur.execute(magic_query, [seq_id])

                # Grab the mRNAs on this sequence
                db_cur.execute('SELECT * FROM tmp_cur_seq WHERE type="mRNA"')
                rnas = db_cur.fetchall()

                for rna in rnas:
                    # Grab the trinotate entry
                    db_cur.execute('SELECT * FROM trinotate WHERE prot_id=? LIMIT 1', [rna[5]])
                    trinotate = db_cur.fetchone()

                    ##### EXON
                    # Grab all exons under this mRNA
                    # TODO grab annotation also
                    db_cur.execute('SELECT start, stop FROM tmp_cur_seq WHERE type="exon" AND parent=?', [rna[0]])
                    rows = db_cur.fetchall()

                    if len(rows) > 0:
                        f.write(str(rows[0][0])+'\t'+str(rows[0][1])+'\texon\n')
                        rows = rows[1:]
                        for row in rows:
                            f.write(str(row[0])+'\t'+str(row[1])+'\n')

                    ## Write the annotations

                    ##### CDS
                    # Grab all CDSs under this mRNA
                    # TODO grab annotation also
                    db_cur.execute('SELECT start, stop FROM tmp_cur_seq WHERE type="CDS" AND parent=?', [rna[0]])
                    rows = db_cur.fetchall()                    

                    if len(rows) > 0:
                        f.write(str(rows[0][0])+'\t'+str(rows[0][1])+'\tCDS\n')
                        rows = rows[1:]
                        for row in rows:
                            f.write(str(row[0])+'\t'+str(row[1])+'\n')
                        cds_ann = parse_gene_ontology(trinotate[8])
                        for annot in cds_ann.split(','):
                            key_val = annot.split('=')
                            if len(key_val) == 2:
                                f.write(key_val[0]+'\t'+key_val[1]+'\n')
                
                db_cur.execute('DROP TABLE tmp_cur_seq')

#!/usr/bin/env python

import sqlite3
import re

def regexp(expr, item):
    reg = re.compile(expr)
    return reg.search(item) is not None

# Takes blast hit and returns it as a CSV of key=value
def parse_blast_hit_gene(blast_hit):
    if len(blast_hit) == 0 or blast_hit == '.':
        return ''

    attributes = ''
    parts = blast_hit.split('`')

    part1_split = parts[0].split('|')
    attributes += "gene="+re.sub(r'_(.*)$', '', part1_split[2])

    return attributes

# Takes blast hit and returns it as a CSV of key=value
def parse_blast_hit_cds(blast_hit):
    if len(blast_hit) == 0 or blast_hit == '.':
        return 'product=hypothetical protein'

    attributes = ''
    parts = blast_hit.split('`')

    matchObj = re.search( r'RecName: Full=([^;]*)', blast_hit)
    if matchObj:
        attributes += "product="+matchObj.group(1)
    else:
        attributes += "product=hypothetical protein"

    return attributes

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
            # TODO try-catch here?
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
        db_conn.create_function('regexp', 2, regexp)

        db_cur = db_conn.cursor()

        with open(fileName, 'w') as f:
            ## Stupid first line??
            f.write('>Feature SeqId\n')

            # Get the sequence ids
            db_cur.execute('SELECT * FROM fasta')

            i = 0
            for seq in db_cur.fetchall():
                i += 1
                if i%100 == 0:
                    print(i)

                f.write('>Feature '+seq[0]+'\n')
                f.write('1\t'+str(len(seq[1]))+'\tREFERENCE\n')
                f.write('\t\t\tPBARC\t12345\n')

                # Create the temporary table
                # db_cur.execute('CREATE TABLE tmp_cur_seq(id TEXT PRIMARY KEY, seq_id TEXT, type TEXT, start INT, stop INT, name TEXT, parent TEXT)')
                
                # TODO: Match columns up with gff columns (i.e. strand isn't last in gff table)
                magic_query = "CREATE TEMP TABLE tmp_cur_seq AS SELECT id, seq_id, type, start, stop, name, parent, strand FROM gff WHERE seq_id=?"

                # Store temporary table of all info pertaining to this sequence
                db_cur.execute(magic_query, [seq[0]])

                # Grab the genes on this sequence
                db_cur.execute('SELECT * FROM tmp_cur_seq WHERE type="gene"')
                genes = db_cur.fetchall()

                for gene in genes:
                    # Find a trinotate entry to give us the gene id for this gene
                    db_cur.execute('SELECT top_blast_hit FROM trinotate WHERE prot_id REGEXP ? LIMIT 1', [gene[5]])
                    trinotate = db_cur.fetchone()

                    # Write the gene
                    f.write(str(gene[3])+'\t'+str(gene[4])+'\tgene\n')

                    gene_ann = 'locus_tag='+gene[5]
                    if trinotate:
                        gene_ann += ','+parse_blast_hit_gene(trinotate[0])
                    for annot in gene_ann.split(','):
                        key_val = annot.split('=')
                        if len(key_val) == 2:
                            f.write('\t\t\t'+key_val[0]+'\t'+key_val[1]+'\n')

                    # Grab the mRNAs on this sequence
                    db_cur.execute('SELECT * FROM tmp_cur_seq WHERE type="mRNA" AND parent=?', [gene[0]])
                    rnas = db_cur.fetchall()

                    for rna in rnas:
                        # Grab the trinotate entry
                        trinotate = None
                        if re.match( r'maker', rna[5]):
                            db_cur.execute('SELECT * FROM trinotate WHERE prot_id=? LIMIT 1', [rna[5]])
                        else:
                            stripped_prot_id = re.sub(r'\.(\d*)', '', rna[5])
                            db_cur.execute('SELECT * FROM trinotate WHERE prot_id=? LIMIT 1', [stripped_prot_id])
                            trinotate = db_cur.fetchone()

                        if not trinotate:
                            print("Couldn't find annotation for "+rna[5]+". Skipping...")
                            continue

                        # Get whether the sequence has start and stop codons
                        db_cur.execute('SELECT EXISTS(SELECT type, name FROM tmp_cur_seq WHERE type="start_codon" AND parent=? LIMIT 1)', [rna[0]])
                        has_start = db_cur.fetchone()[0]
                        db_cur.execute('SELECT EXISTS(SELECT type, name FROM tmp_cur_seq WHERE type="stop_codon" AND parent=? LIMIT 1)', [rna[0]])
                        has_stop = db_cur.fetchone()[0]

                        # Annotation stuff
                        split_prot_id = rna[5].split('_')
                        cds_exon_ann_base = ''
                        cds_exon_ann_base += 'protein_id='+'gnl|PBARC|'+rna[5]+','
                        cds_exon_ann_base += 'transcript_id='+'gnl|PBARC|'+split_prot_id[0]+'_mrna'+split_prot_id[1]+','
                        cds_exon_ann_base += parse_blast_hit_cds(trinotate[3])

                        ##### EXON
                        # Grab all exons under this mRNA
                        db_cur.execute('SELECT start, stop, strand FROM tmp_cur_seq WHERE type="exon" AND parent=?', [rna[0]])
                        rows = db_cur.fetchall()

                        ######### TODO: HEY BRIAN I'M WRITING COORDINATES
                        if len(rows) > 0:
                            strand = rows[0][2]

                            # Sort dem lists
                            lrows = list()
                            for row in rows:
                                lrows.append(list(row[:-1]))
                            rows = lrows

                            rows.sort()
                            for row in rows:
                                row.sort()

                            if strand == '-':
                                rows.reverse()
                                for row in rows:
                                    row.reverse()

                            if has_start == 0:
                                f.write('<')
                            if has_stop == 0 and len(rows) == 1:
                                f.write(str(rows[0][0])+'\t>'+str(rows[0][1])+'\tmRNA\n') # Exons are actually mRNA. What?
                            else:
                                f.write(str(rows[0][0])+'\t'+str(rows[0][1])+'\tmRNA\n') # Exons are actually mRNA. What?
                            rows = rows[1:]
                            for row in rows[:-1]:
                                f.write(str(row[0])+'\t'+str(row[1])+'\n')
                            if len(rows) > 0:
                                if has_stop == 0:
                                    f.write(str(rows[len(rows)-1][0])+'\t>'+str(rows[len(rows)-1][1])+'\n')
                                else:
                                    f.write(str(rows[len(rows)-1][0])+'\t'+str(rows[len(rows)-1][1])+'\n')

                            exon_ann = cds_exon_ann_base
                            for annot in exon_ann.split(','):
                                key_val = annot.split('=')
                                if len(key_val) == 2:
                                    f.write('\t\t\t'+key_val[0]+'\t'+key_val[1]+'\n')

                        ## Write the annotations

                        ##### CDS
                        # Grab all CDSs under this mRNA
                        db_cur.execute('SELECT start, stop, strand FROM tmp_cur_seq WHERE type="CDS" AND parent=?', [rna[0]])
                        rows = db_cur.fetchall()                    

                        ######### TODO: HEY BRIAN I'M WRITING COORDINATES
                        if len(rows) > 0:
                            strand = rows[0][2]

                            # Sort dem lists
                            lrows = list()
                            for row in rows:
                                lrows.append(list(row[:-1]))
                            rows = lrows

                            rows.sort()
                            for row in rows:
                                row.sort()

                            if strand == '-':
                                rows.reverse()
                                for row in rows:
                                    row.reverse()

                            if has_start == 0:
                                f.write('<')
                            if has_stop == 0 and len(rows) == 1:
                                f.write(str(rows[0][0])+'\t>'+str(rows[0][1])+'\tCDS\n')
                            else:
                                f.write(str(rows[0][0])+'\t'+str(rows[0][1])+'\tCDS\n')
                            rows = rows[1:]
                            for row in rows[:-1]:
                                f.write(str(row[0])+'\t'+str(row[1])+'\n')
                            if len(rows) > 0:
                                if has_stop == 0:
                                    f.write(str(rows[len(rows)-1][0])+'\t>'+str(rows[len(rows)-1][1])+'\n')
                                else:
                                    f.write(str(rows[len(rows)-1][0])+'\t'+str(rows[len(rows)-1][1])+'\n')

                            cds_ann = cds_exon_ann_base+','+parse_gene_ontology(trinotate[8])
                            for annot in cds_ann.split(','):
                                key_val = annot.split('=')
                                if len(key_val) == 2:
                                    f.write('\t\t\t'+key_val[0]+'\t'+key_val[1]+'\n')
                
                db_cur.execute('DROP TABLE tmp_cur_seq')

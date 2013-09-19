#!/usr/bin/env python

import sqlite3
import re
from feature_tbl_entry import FeatureTblEntry

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

	element = element.replace(',', '')
        vals = element.split('^')
        if vals[1] == 'molecular_function':
            attributes += 'go_function='
        elif vals[1] == 'cellular_component':
            attributes += 'go_component='
        elif vals[1] == 'biological_process':
            attributes += 'go_process='
        attributes += vals[2]+'|'+vals[0][3:]+'||IEA'
    return attributes

def reverse_indices(row):
    tmp = row[0]
    row[0] = row[1]
    row[1] = tmp
    return row


class FeatureTblWriter:

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
                magic_query = "CREATE TEMP TABLE tmp_cur_seq AS SELECT id, seq_id, type, start, stop, name, parent, strand, phase FROM gff WHERE seq_id=?"

                # Store temporary table of all info pertaining to this sequence
                db_cur.execute(magic_query, [seq[0]])

                # Grab the genes on this sequence
                db_cur.execute('SELECT * FROM tmp_cur_seq WHERE type="gene"')
                genes = db_cur.fetchall()

                for gene in genes:
                    # Whether or not the gene is good to write
                    gene_good = True
                    exon_entries = []
                    cds_entries = []

                    # Find a trinotate entry to give us the gene id for this gene
                    db_cur.execute('SELECT top_blast_hit FROM trinotate WHERE prot_id REGEXP ? LIMIT 1', [gene[5]])
                    trinotate = db_cur.fetchone()

                    # Write the gene
                    gene_entry = FeatureTblEntry()
                    gene_entry.set_type('gene')
                    gene_entry.set_strand(gene[7])
                    gene_entry.set_partial_info(True, True)
                    if gene[8] != '.':
                        gene_entry.set_phase(int(gene[8]))
                    gene_entry.add_coordinates(gene[3], gene[4])
                    gene_ann = 'locus_tag='+gene[5]
                    if trinotate:
                        gene_ann += ','+parse_blast_hit_gene(trinotate[0])
                    for annot in gene_ann.split(','):
                        key_val = annot.split('=')
                        if len(key_val) == 2:
                            gene_entry.add_annotation(key_val[0], key_val[1])
                    if gene_entry.get_total_length() < 150:
                        gene_good = False

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
                        db_cur.execute('SELECT start, stop, strand, phase FROM tmp_cur_seq WHERE type="exon" AND parent=?', [rna[0]])
                        exons = db_cur.fetchall()

                        if len(exons) > 0:
                            entry = FeatureTblEntry()
                            entry.set_type('mRNA') # Exons are actually mRNA?
                            entry.set_partial_info(has_start, has_stop)
                            entry.set_strand(exons[0][2])
                            if exons[0][3] != '.':
                                entry.set_phase(int(exons[0][3]))
                            [entry.add_coordinates(exon[0], exon[1]) for exon in exons]
                            exon_ann = cds_exon_ann_base
                            for annot in exon_ann.split(','):
                                key_val = annot.split('=')
                                if len(key_val) == 2:
                                    entry.add_annotation(key_val[0], key_val[1])
                            exon_entries.append(entry)
                            if entry.get_total_length() < 150:
                                gene_good = False
                        else:
                            gene_good = False

                        ## Write the annotations

                        ##### CDS
                        # Grab all CDSs under this mRNA
                        db_cur.execute('SELECT start, stop, strand, phase FROM tmp_cur_seq WHERE type="CDS" AND parent=?', [rna[0]])
                        CDSs = db_cur.fetchall()                        

                        if len(CDSs) > 0:
                            entry = FeatureTblEntry()
                            entry.set_type('CDS')
                            entry.set_partial_info(has_start, has_stop)
                            entry.set_strand(CDSs[0][2])
                            if CDSs[0][3] != '.':
                                entry.set_phase(int(CDSs[0][3]))
                            [entry.add_coordinates(CDS[0], CDS[1]) for CDS in CDSs]
                            cds_ann = cds_exon_ann_base+','+parse_gene_ontology(trinotate[8])
                            for annot in cds_ann.split(','):
                                key_val = annot.split('=')
                                if len(key_val) == 2:
                                    entry.add_annotation(key_val[0], key_val[1])
                            cds_entries.append(entry)
                            if entry.get_total_length() < 150:
                                gene_good = False
                        else:
                            gene_good = False

                    if len(exon_entries) != len(cds_entries):
                        gene_good = False
                    if gene_good:
                        f.write(gene_entry.write_to_string())
                        for i in range(len(exon_entries)):
                            f.write(exon_entries[i].write_to_string())
                            f.write(cds_entries[i].write_to_string())
                    else:
                        print("Removed bad gene ", gene[5])

                db_cur.execute('DROP TABLE tmp_cur_seq')

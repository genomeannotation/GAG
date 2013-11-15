#!/usr/bin/env python

import re
from feature_tbl_entry import FeatureTblEntry

# Takes blast hit and returns it as a list of key value pairs
def parse_blast_hit_gene(blast_hit):
    if len(blast_hit) == 0 or blast_hit == '.':
        return ''

    attributes = []
    parts = blast_hit.split('`')

    part1_split = parts[0].split('|')
    attributes.append(["gene", re.sub(r'_(.*)$', '', part1_split[2])])

    return attributes

# Takes blast hit and returns it as a list of key value pairs
def parse_blast_hit_cds(blast_hit):
    if len(blast_hit) == 0 or blast_hit == '.':
        return [['product', 'hypothetical protein']]

    attributes = []
    parts = blast_hit.split('`')

    matchObj = re.search( r'RecName: Full=([^;]*)', blast_hit)
    if matchObj:
        attributes.append(["product", matchObj.group(1)])
    else:
        attributes.append(["product", "hypothetical protein"])

    return attributes

# Takes gene ontology and returns it as a list of key value pairs
def parse_gene_ontology(gene_ontology):
    if len(gene_ontology) == 0 or gene_ontology == '.':
        return ''

    attributes = []
    elements = gene_ontology.split('`')
    for element in elements:
        element = element.replace(',', '')
        vals = element.split('^')
        key = ''
        if vals[1] == 'molecular_function':
            key = 'go_function'
        elif vals[1] == 'cellular_component':
            key = 'go_component'
        elif vals[1] == 'biological_process':
            key = 'go_process'
        attributes.append([key, vals[2]+'|'+vals[0][3:]+'||IEA'])
    return attributes

class Annotator:

    def __init__(self):
        self.entries = []

    def read_from_file(self, filename):
        with open(filename, "r") as f:
            for line in f:
                if line[0] == '#':
                    continue

                self.entries.append(line.split('\t'))

    def annotate_gene(self, gene):
        gene.add_annotation('locus_tag', gene.name)
        for entry in self.entries:
            if entry[2][:-3] == gene.name: # Chop off the -RA stuff
                gene.add_annotations(parse_blast_hit_gene(entry[3]))
                return

    def annotate_cds(self, cds):
        split_prot_id = cds.name.split('_')
        cds.add_annotation('protein_id', 'gnl|PBARC|'+cds.name)
        cds.add_annotation('transcript_id', 'gnl|PBARC|'+split_prot_id[0]+'_mrna'+split_prot_id[1])

        for entry in self.entries:
            if entry[2] == cds.name:
                cds.add_annotations(parse_blast_hit_cds(entry[3]))
                cds.add_annotations(parse_gene_ontology(entry[8]))
                return

        cds.add_annotation('product', 'hypothetical protein')

    def annotate_mrna(self, mrna):
        split_prot_id = mrna.name.split('_')                
        mrna.add_annotation('protein_id', 'gnl|PBARC|'+mrna.name)
        mrna.add_annotation('transcript_id', 'gnl|PBARC|'+split_prot_id[0]+'_mrna'+split_prot_id[1])

        for entry in self.entries:
            if entry[2] == mrna.name:
                mrna.add_annotations(parse_blast_hit_cds(entry[3]))
                return

        mrna.add_annotation('product', 'hypothetical protein')


#!/usr/bin/env python

import sys
from gene_part import GenePart, CDS, Exon
from mrna import MRNA
from gene import Gene

# TODO apply_bed(bed) -- calls gene.trim(bed.get_coords) for all genes
# where bed.contains(gene.seq_id)
# damn i just about wrote it already...

class GFF:

    def __init__(self):
        self.genes = []
        self.current_gene = None
        self.current_mrna = None
        self.current_exon = None
        self.current_cds = None

    def __str__(self):
        result = "GFF containing "
        result += str(len(self.genes))
        result += " genes\n"
        return result

    def validate_line(self, line):
        if len(line) is not 9:
            return False
        return True

    def line_type(self, line):
        return line[2]

    def validate_first_line(self, line):
        if self.line_type(line) is 'gene':
            return self.validate_line(line)
        else:
            return False

    # returns a dict with id, name, parent_id (if present)
    def parse_attributes(self, attr):
        result = {}
        split_attr = attr.split(';')
        # make sure that worked :)
        if len(split_attr) < 2:
            sys.stderr.write('error trying to parse "' + attr + '"\n')
            return None
        try:
            result['identifier'] = split_attr[0].split('=')[1]
            result['name'] = split_attr[1].split('=')[1]
        except IndexError:
            sys.stderr.write('error trying to parse "' + attr + '"\n')
            return None
        if len(split_attr) is 3:
            try:
                result['parent_id'] = split_attr[2].split('=')[1]
            except IndexError:
                sys.stderr.write('error trying to parse "' + attr + '"\n')
                if self.current_gene:
                    sys.stderr.write('occurred at ' + str(self.current_gene))
                return None
        return result
        

    def extract_cds_args(self, line):
        result = {'indices': [int(line[3]), int(line[4])], 'phase': int(line[7])}
        if isinstance(line[7], float):
            result['score'] = line[7]
        attribs = self.parse_attributes(line[8])
        result.update(attribs)
        return result

    def extract_exon_args(self, line):
        result = {'indices': [int(line[3]), int(line[4])], 'score': line[5]}
        attribs = self.parse_attributes(line[8])
        result.update(attribs)
        return result

    def extract_other_feature_args(self, line):
        result = {'feature_type': line[2], 'indices': [int(line[3]), int(line[4])]}
        attribs = self.parse_attributes(line[8])
        result.update(attribs)
        return result

    def extract_mrna_args(self, line):
        result = {'indices': [int(line[3]), int(line[4])]}
        attribs = self.parse_attributes(line[8])
        result.update(attribs)
        return result        

    def extract_gene_args(self, line):  # TODO 'score'
        result = {'seq_name': line[0], 'source': line[1], 'indices': [int(line[3]), int(line[4])], 'strand': line[6]}
        attribs = self.parse_attributes(line[8])
        result.update(attribs)
        return result

    def remove_first_cds_segment_if_shorter_than(self, min_length):
        if self.genes:
            for gene in self.genes:
                gene.remove_first_cds_segment_if_shorter_than(min_length)

    def remove_mrnas_with_cds_shorter_than(self, min_length):
        if self.genes:
            to_remove = []
            for gene in self.genes:
                gene.remove_mrnas_with_cds_shorter_than(min_length)
                if not gene.mrnas:
                    to_remove.append(gene)
            for g in to_remove:
                self.genes.remove(g)

    def update_cds(self, line):
        args = self.extract_cds_args(line)
        self.current_cds.add_indices(args['indices'])
        self.current_cds.add_phase(args['phase'])
        self.current_cds.add_identifier(args['identifier'])
        self.current_cds.add_name(args['name'])
        if 'score' in args:
            self.current_cds.add_score(args['score'])

    def update_exon(self, line):
        args = self.extract_exon_args(line)
        self.current_exon.add_indices(args['indices'])
        self.current_exon.add_identifier(args['identifier'])
        self.current_exon.add_name(args['name'])
        if 'score' in args:
            self.current_exon.add_score(args['score'])

    def process_line(self, line):
        ltype = self.line_type(line)
        if ltype == 'gene':
            self.process_gene_line(line)
        elif ltype == 'mRNA':
            self.process_mrna_line(line)
        elif ltype == 'CDS':
            self.process_cds_line(line)
        elif ltype == 'exon':
            self.process_exon_line(line)
        else:
            self.process_other_feature_line(line)

    def process_gene_line(self, line):
        if self.current_gene:
            self.wrap_up_gene()
            self.process_gene_line(line)
        else:
            kwargs = self.extract_gene_args(line)
            self.current_gene = Gene(**kwargs)

    def process_mrna_line(self, line):
        if self.current_mrna:
            self.wrap_up_mrna()
            self.process_mrna_line(line)
        else:
            kwargs = self.extract_mrna_args(line)
            self.current_mrna = MRNA(**kwargs)

    def process_cds_line(self, line):
        if self.current_cds:
            self.update_cds(line)
        else:
            kwargs = self.extract_cds_args(line)
            self.current_cds = CDS(**kwargs)

    def process_exon_line(self, line):
        if self.current_exon:
            self.update_exon(line)
        else:
            kwargs = self.extract_exon_args(line)
            self.current_exon = Exon(**kwargs)

    def process_other_feature_line(self, line):
        pass

    def wrap_up_gene(self):
        if self.current_mrna:
            self.wrap_up_mrna()
        self.genes.append(self.current_gene)
        self.current_gene = None

    def wrap_up_mrna(self):
        if self.current_cds:
            # TODO check parent_id?
            self.current_mrna.set_cds(self.current_cds)
            self.current_cds = None
        if self.current_exon:
            # TODO check parent_id? may have to change tests...
            self.current_mrna.set_exon(self.current_exon)
            self.current_exon = None
        self.current_gene.add_mrna(self.current_mrna)
        self.current_mrna = None

    def read_file(self, reader):
        for line in reader:
            if len(line) == 0 or line[0].startswith('#'):
                continue
            else:
                if self.validate_line(line):
                    self.process_line(line)
        self.wrap_up_gene()

    def apply_bed(self, bed):
        for gene in self.genes:
            if bed.contains(gene.seq_name):
                coords = bed.get_coordinates(gene.seq_name)
                gene.trim(coords)

    def subset_gff(self, seqlist):
        self.genes = [g for g in self.genes if g.seq_name in seqlist]

    def remove_empty_genes(self):
        self.genes = [g for g in self.genes if not g.is_empty()]

    def remove_all_gene_segments(self, prefix):
        if len(prefix) > 0:
            self.genes = [g for g in self.genes if not prefix in g.name]

    # takes list of mrna names, returns list of gene names
    def get_mrnas_parent_gene_names(self, mrnalist):
        gene_names = []
        for gene in self.genes:
            for mrna in mrnalist:
                if gene.contains_mrna_named(mrna) and gene.name not in gene_names:
                    gene_names.append(gene.name)
        return gene_names

    def gene_name_to_prefix(self, gene_name):
        return gene_name.split('.')[0]

    def prefix_match(self, gene, prefixes):
        for prefix in prefixes:
            if prefix in gene.name:
                return True
        return False

    def remove_genes_by_prefixes(self, prefixes):
        self.genes = [g for g in self.genes if not self.prefix_match(g, prefixes)]


    def obliterate_genes_related_to_mrnas(self, mrna_names):
        parent_genes = self.get_mrnas_parent_gene_names(mrna_names)
        prefixes = [self.gene_name_to_prefix(g) for g in parent_genes]
        self.remove_genes_by_prefixes(prefixes)

    def invalidate_region(self, seq, start, stop):
        for gene in self.genes:
            if gene.seq_name == seq:
                gene.invalidate_region(start, stop)

    def contains_gene_on_seq(self, seq_id):
        for gene in self.genes:
            if gene.seq_name == seq_id:
                return True
        return False

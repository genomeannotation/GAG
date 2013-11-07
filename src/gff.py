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
        result += " genes"
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
        if not self.current_mrna:
            sys.stderr.write("trying to add a feature but no mRNA to add it to... here is the line in question: " + str(line))
        else:
            kwargs = self.extract_other_feature_args(line)
            feat = GenePart(**kwargs)
            self.current_mrna.other_features.append(feat)

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
            if bed.contains(gene.seq_id):
                coords = bed.get_coordinates(gene.seq_id)
                gene.trim(coords)
                    





#!/usr/bin/env python

import sys
from feature_classes import GenePart, CDS, Exon, MRNA, Gene

class GFF:

    def __init__(self):
        self.genes = []
        self.current_gene = None
        self.current_mrna = None
        self.current_exon = None
        self.current_cds = None

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
                return None
        return result
        

    def extract_cds_args(self, line):
        result = {'indices': [int(line[3]), int(line[4])], 'phase': int(line[7])}
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
        

    # returns a dict 
    def extract_gene_args(self, line):  # TODO 'score'
        result = {'seq_name': line[0], 'source': line[1], 'indices': [int(line[3]), int(line[4])], 'strand': line[6]}
        attribs = self.parse_attributes(line[8])
        result.update(attribs)
        return result

    def update_cds(self, line):
        # TODO
        print "foo"

    def process_cds_line(self, line):
        if self.current_cds:
            self.update_cds(line)
        else:
            args = self.extract_cds_args(line)





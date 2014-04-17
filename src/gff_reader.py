#!/usr/bin/env python

import sys
import copy
from src.gene_part import GenePart, CDS, Exon
from src.mrna import MRNA
from src.gene import Gene

class GFFReader:

    def __init__(self):
        self.genes = {}
        self.mrnas = {}
        self.orphans = []
        self.give_up = False # we are strong for now
        self.skipped_features = 0

    def validate_line(self, line):
        """Returns list of fields if valid, empty list if not."""
        splitline = line.split('\t')
        if len(splitline) is not 9:
            return []
        if not "ID" in splitline[8]:
            return []
        # Everything except genes must have parent id
        if not "Parent" in splitline[8] and not splitline[2] == "gene":
            return []
        return splitline

    def line_type(self, line):
        """Returns type of feature, as denoted by 3rd field in list."""
        return line[2]

    def parse_attributes(self, attr):
        """Returns a dict with id and parent_id (if present)
        if not, returns empty dict
        Also adds annotations if present
        """
        result = {}
        annotations = []
        # Sanitize and split attributes up
        split_attr = attr.strip(' \t\n;').split(';')
        for pair in split_attr:
            splitpair = pair.split('=')
            if len(splitpair) != 2:
                continue
            if splitpair[0] == "ID":
                result['identifier'] = splitpair[1]
            elif splitpair[0] == "Parent":
                result['parent_id'] = splitpair[1]
            elif splitpair[0] == "Dbxref" or splitpair[0] == "Ontology_term":
                annotations.append(splitpair)
        # Make sure we found an ID
        if "identifier" not in result:
            return {}
        # Add annotations if we found any
        if annotations:
            result["annotations"] = annotations
        return result

    def extract_cds_args(self, line):
        result = {'indices': [int(line[3]), int(line[4])], \
                'strand': line[6], 'phase': int(line[7])}
        if isinstance(line[7], float):
            result['score'] = line[7]
        attribs = self.parse_attributes(line[8])
        
        if not attribs:
            return None

        result.update(attribs)
        return result

    def extract_exon_args(self, line):
        result = {'indices': [int(line[3]), int(line[4])], 'strand': line[6]}
        if line[5] != '.':
            result['score'] = float(line[5])
        attribs = self.parse_attributes(line[8])

        if not attribs:
            return None

        result.update(attribs)
        return result

    def extract_mrna_args(self, line):
        result = {'indices': [int(line[3]), int(line[4])]}
        attribs = self.parse_attributes(line[8])

        if not attribs:
            return None

        result.update(attribs)
        return result        

    def extract_gene_args(self, line):  
        result = {'seq_name': line[0], 'source': line[1], \
                  'indices': [int(line[3]), int(line[4])], 'strand': line[6]}
        attribs = self.parse_attributes(line[8])

        if not attribs:
            return None

        result.update(attribs)
        return result

    def extract_other_feature_args(self, line):
        result = {'feature_type': line[2], 'indices': [int(line[3]), int(line[4])]}
        attribs = self.parse_attributes(line[8])
        result.update(attribs)
        return result

    def update_cds(self, line, cds):
        args = self.extract_cds_args(line)
        cds.add_indices(args['indices'])
        cds.add_phase(args['phase'])
        cds.add_identifier(args['identifier'])
        if 'score' in args:
            cds.add_score(args['score'])

    def update_exon(self, line, exon):
        args = self.extract_exon_args(line)
        exon.add_indices(args['indices'])
        exon.add_identifier(args['identifier'])
        if 'score' in args:
            exon.add_score(args['score'])

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
        elif ltype == 'start_codon' or ltype == 'stop_codon':
            self.process_other_feature_line(line)
	else:
            self.skipped_features += 1

    def process_gene_line(self, line):
        kwargs = self.extract_gene_args(line)
        if not kwargs:
            return
        gene_id = kwargs['identifier']
        self.genes[gene_id] = Gene(**kwargs)

    def process_mrna_line(self, line):
        kwargs = self.extract_mrna_args(line)
        if not kwargs:
            return
        mrna_id = kwargs['identifier']
        self.mrnas[mrna_id] = MRNA(**kwargs)

    def process_cds_line(self, line):
        kwargs = self.extract_cds_args(line)
        if not kwargs:
            return
        parent_id = kwargs['parent_id']
        if parent_id not in self.mrnas:
            self.orphans.append(line)
            return
        parent_mrna = self.mrnas[parent_id]
        if parent_mrna.cds:
            self.update_cds(line, parent_mrna.cds)
        else:
            parent_mrna.cds = CDS(**kwargs)

    def process_exon_line(self, line):
        kwargs = self.extract_exon_args(line)
        if not kwargs:
            return
        parent_id = kwargs['parent_id']
        if parent_id not in self.mrnas:
            self.orphans.append(line)
            return
        parent_mrna = self.mrnas[parent_id]
        if parent_mrna.exon:
            self.update_exon(line, parent_mrna.exon)
        else:
            parent_mrna.exon = Exon(**kwargs)

    def process_other_feature_line(self, line):
        kwargs = self.extract_other_feature_args(line)
        if not kwargs:
            return
        parent_id = kwargs['parent_id']
        if parent_id not in self.mrnas:
            self.orphans.append(line)
            return
        parent_mrna = self.mrnas[parent_id]
        parent_mrna.other_features.append(GenePart(**kwargs))

    def read_file(self, reader):
        # First pass, pulling out all genes and mRNAs
        #  and placing child features if possible
        for line in reader:
            if len(line) == 0 or line[0].startswith('#'):
                continue
            if self.give_up:
                return
            splitline = self.validate_line(line)
            if splitline:
                self.process_line(splitline)

        # Second pass, placing child features which 
        # preceded their parents in the first pass
        orphans = copy.deepcopy(self.orphans)
        for splitline in orphans:
            self.process_line(splitline)
           
        # Add mRNAs to their parent genes
        for mrna in self.mrnas.values():
            parent_gene = self.genes[mrna.parent_id]
            parent_gene.mrnas.append(mrna)

        if self.skipped_features > 0:
            sys.stderr.write("Warning: skipped "+str(self.skipped_features)+" uninteresting features.")
        return self.genes.values()


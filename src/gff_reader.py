#!/usr/bin/env python

import sys
import traceback
from src.gene_part import CDS, Exon
from src.mrna import MRNA
from src.gene import Gene

class GFFReader:

    def __init__(self):
        self.genes = []
        self.current_gene = None
        self.current_mrna = None
        self.current_exon = None
        self.current_cds = None
        self.current_line = 0 # Not even reading a file yet
        self.give_up = False # we are strong for now

    def validate_line(self, line):
        """Returns list of fields if valid, empty list if not."""
        splitline = line.split('\t')
        if len(splitline) is not 9:
            return []
        if not 'ID' in splitline[8]:
            return []
        return splitline

    def line_type(self, line):
        """Returns type of feature, as denoted by 3rd field in list."""
        return line[2]

    # returns a dict with id and parent_id (if present)
    # if not, returns empty dict
    def parse_attributes(self, attr):
        split_attr = attr.split(';')
        try:
            keys = [val.split('=')[0] for val in split_attr]
            vals = [val.split('=')[1] for val in split_attr]
        except IndexError as ie:
            sys.stderr.write("IndexError trying to split attributes: " + str(split_attr))
            return None

        
        attr_dict = dict(zip(keys, vals)) # Our parameter dictionary
        
        result = {}

        try:
            result['identifier'] = attr_dict['ID'].strip()
            if 'Parent' in attr_dict:
                result['parent_id'] = attr_dict['Parent'].strip()
        except KeyError as ke:
            print("\nError reading GFF mRNA entry at line "+str(self.current_line)+": required attribute '"+ke.args[0]+"' doesn't exist.\n")

            go_on = raw_input("\n\nAttempt to continue? (y/n): ")
            if go_on != 'y' and go_on != 'Y': # Didn't select Y, get outta here!
                self.give_up = True

            return {}

        return result
        

    def extract_cds_args(self, line):
        result = {'indices': [int(line[3]), int(line[4])], \
                  'phase': int(line[7])}
        if isinstance(line[7], float):
            result['score'] = line[7]
        attribs = self.parse_attributes(line[8])
        
        if not attribs:
            return None

        result.update(attribs)
        return result

    def extract_exon_args(self, line):
        result = {'indices': [int(line[3]), int(line[4])], 'score': float(line[5])}
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
        if 'score' in args:
            self.current_cds.add_score(args['score'])

    def update_exon(self, line):
        args = self.extract_exon_args(line)
        self.current_exon.add_indices(args['indices'])
        self.current_exon.add_identifier(args['identifier'])
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

            if not kwargs:
                return

            self.current_gene = Gene(**kwargs)

    def process_mrna_line(self, line):
        if self.current_mrna:
            self.wrap_up_mrna()
            self.process_mrna_line(line)
        else:
            kwargs = self.extract_mrna_args(line)

            if not kwargs:
                return

            self.current_mrna = MRNA(**kwargs)

    def process_cds_line(self, line):
        if self.current_cds:
            self.update_cds(line)
        else:
            kwargs = self.extract_cds_args(line)

            if not kwargs:
                return

            self.current_cds = CDS(**kwargs)

    def process_exon_line(self, line):
        if self.current_exon:
            self.update_exon(line)
        else:
            kwargs = self.extract_exon_args(line)

            if not kwargs:
                return

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
            self.current_mrna.set_cds(self.current_cds)
            self.current_cds = None
        if self.current_exon:
            self.current_mrna.set_exon(self.current_exon)
            self.current_exon = None
        self.current_gene.add_mrna(self.current_mrna)
        self.current_mrna = None

    def read_file(self, reader):
        self.current_line = 0 # aaaand begin!
        for line in reader:
            if len(line) == 0:
                continue
            if self.give_up:
                return

            self.current_line += 1

            try:
                if len(line) == 0 or line[0].startswith('#'):
                    continue
                else:
                    if self.validate_line(line):
                        self.process_line(line)
            except:
                print("\nException raised while reading GFF line: "+str(self.current_line)+"\n\n")
                print("The line looks like this:\n")
                print(line)
                print(traceback.format_exc())
                go_on = raw_input("\n\nAttempt to continue? (y/n): ")
                if go_on != 'y' and go_on != 'Y': # Didn't select Y, get outta here!
                    return
        self.wrap_up_gene()

    def subset_gff(self, seqlist):
        self.genes = [g for g in self.genes if g.seq_name in seqlist]

    def remove_empty_genes(self):
        self.genes = [g for g in self.genes if not g.is_empty()]

    def remove_all_gene_segments(self, prefix):
        if len(prefix) > 0:
            self.genes = [g for g in self.genes if not prefix in g.identifier]

    def prefix_match(self, gene, prefixes):
        for prefix in prefixes:
            if prefix in gene.identifier:
                return True
        return False

    def remove_genes_by_prefixes(self, prefixes):
        self.genes = \
                [g for g in self.genes if not self.prefix_match(g, prefixes)]


    def remove_genes_marked_for_removal(self):
        for gene in reversed(self.genes):
            if gene.indices[0] == 0 and gene.indices[1] == 0:
                self.genes.remove(gene)

    def invalidate_region(self, seq, start, stop):
        for gene in self.genes:
            if gene.seq_name == seq:
                gene.invalidate_region(start, stop)

    def contains_gene_on_seq(self, seq_id):
        for gene in self.genes:
            if gene.seq_name == seq_id:
                return True
        return False

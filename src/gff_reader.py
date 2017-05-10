#!/usr/bin/env python
# coding=utf-8

import sys
import copy
from src.gene_part import GenePart
from src.cds import CDS
from src.exon import Exon
from src.xrna import XRNA
from src.gene import Gene


def has_multiple_parents(attr):
    split_attr = attr.split(";")
    for attr in split_attr:
        if "Parent" in attr:
            parent_id_field = attr.split("=")[1]
            if "," in parent_id_field:
                return True
    return False


def remove_parent_info_from_attr(fields):
    result = []
    for field in fields:
        if "Parent" not in field:
            result.append(field)
    return result


def get_parents_from_list_of_attributes(fields):
    """Returns a list of parent ids from a list of column 9 entries."""
    for field in fields:
        if "Parent" in field:
            return field.split("=")[1].split(",")


def split_multi_parent_line(fields):
    """Returns a list of lines, one for each parent in a multiparent line."""
    split_attr = fields[8].split(";")
    parents = get_parents_from_list_of_attributes(split_attr)
    attr_without_parents = remove_parent_info_from_attr(split_attr)
    all_lines = []
    # import pdb; pdb.set_trace()
    for parent in parents:
        line = fields[:8]
        new_attr_list = copy.deepcopy(attr_without_parents)
        new_attr_list.append("Parent=" + parent)
        new_attr_string = ";".join(new_attr_list)
        line.append(new_attr_string)
        all_lines.append(line)
    return all_lines


class GFFReader(object):
    def __init__(self):
        self.genes = {}
        self.mrnas = {}
        self.non_genes = []
        self.orphans = []
        self.skipped_features = 0

    @staticmethod
    def validate_line(line):
        """Returns list of lists of fields if valid, empty list if not.

        List of lists of fields -- because lines with multiple parents
        are split into multiple lines."""
        splitline = line.split('\t')
        if len(splitline) is not 9:
            print("not enough columns:" + line)
            return []
        if not "ID" in splitline[8]:
            print("No ID:" + line)
            return []
        if not int(splitline[3]) <= int(splitline[4]):
            print("stop greater than start:" + line)
            return []
        # Special case: mature miRNA and uORF have "Derives_from" relationship
        if splitline[2] in ("uORF", "miRNA") and "Derives_from" in splitline[8]:
            pass
        # Special case: IES and CBS features have no parents
        if splitline[2] in ("chromosome_breakage_sequence", "internal_eliminated_sequence"):
            pass
        # Everything except genes must have parent id
        elif not "Parent" in splitline[8] and \
           not (splitline[2] in ("gene", "pseudogene", "transposable_element_gene")):
            print("no parent:" + line)
            return []
        if has_multiple_parents(splitline[8]):
            splitlines = split_multi_parent_line(splitline)
            return splitlines
        else:
            return [splitline]

    @staticmethod
    def line_type(line):
        """Returns type of feature, as denoted by 3rd field in list."""
        return line[2]

    @staticmethod
    def parse_attributes(attr, ftype=""):
        """Returns a dict with id, name and parent_id (if present)

        If not, returns empty dict
        Also adds annotations if present
        """
        result = {}
        annotations = {}
        # Sanitize and split attributes up
        split_attr = attr.strip(' \t\n;').split(';')
        for pair in split_attr:
            splitpair = pair.split('=')
            if len(splitpair) != 2:
                continue
            key = splitpair[0]
            value = splitpair[1]
            if key in "ID":
                result['identifier'] = value
            elif key == "Name":
                result['name'] = value
            elif key == "Parent":
                result['parent_id'] = value
            elif key == "Derives_from":
                result['derives_from'] = value
            elif (key == "Dbxref" or
                  key == "Ontology_term" or
                  key == "product"):
                if key in annotations.keys():
                    # allow for annotations in the style of "Dbxref=PFAM:foo,PRINTS:bar"
                    annotations[key].extend(value.split(','))
                else:
                    annotations[key] = value.split(',')  # always a list :)
        # Make sure we found an ID
        if "identifier" not in result:
            return {}
        # Add annotations if we found any
        if annotations:
            result["annotations"] = annotations
        # Delete gene name if it is the same as the ID
        if 'name' in result and result['identifier'] == result['name']:
            del result['name']
        # Extract parent_id from derives_from identifier
        if not ftype.endswith("gene") and 'parent_id' not in result and 'derives_from' in result:
            result['parent_id'] = result['derives_from'].split(".")[0]
        return result

    def extract_cds_args(self, line):
        """Pulls CDS arguments from a gff line and returns them in a dictionary."""
        result = {'indices': [int(line[3]), int(line[4])],
                  'strand': line[6], 'phase': int(line[7])}
        if isinstance(line[7], float):
            result['score'] = line[7]
        attribs = self.parse_attributes(line[8], ftype=line[2])

        if not attribs:
            return None

        if "annotations" in attribs:
            del attribs["annotations"]

        result.update(attribs)
        return result

    def extract_exon_args(self, line):
        """Pulls Exon arguments from a gff line and returns them in a dictionary."""
        result = {'indices': [int(line[3]), int(line[4])], 'strand': line[6]}
        if line[5] != '.':
            result['score'] = float(line[5])
        attribs = self.parse_attributes(line[8], ftype=line[2])

        if not attribs:
            return None

        if "annotations" in attribs:
            del attribs["annotations"]

        result.update(attribs)
        return result

    def extract_mrna_args(self, line):
        """Pulls XRNA arguments from a gff line and returns them in a dictionary."""
        result = {'indices': [int(line[3]), int(line[4])], 'strand': line[6],
                  'seq_name': line[0], 'source': line[1]}
        attribs = self.parse_attributes(line[8], ftype=line[2])

        if not attribs:
            return None

        if 'name' in attribs:
            del attribs['name']

        result.update(attribs)
        return result

    def extract_gene_args(self, line):
        """Pulls Gene arguments from a gff line and returns them in a dictionary."""
        result = {'seq_name': line[0], 'source': line[1],
                  'indices': [int(line[3]), int(line[4])], 'strand': line[6]}
        attribs = self.parse_attributes(line[8], ftype=line[2])

        if not attribs:
            return None

        result.update(attribs)
        return result

    def extract_other_feature_args(self, line):
        """Pulls GenePart arguments from a gff line and returns them in a dictionary."""
        result = {'feature_type': line[2], 'indices': [int(line[3]), int(line[4])]}
        attribs = self.parse_attributes(line[8], ftype=line[2])

        if not attribs:
            return None

        if "annotations" in attribs:
            del attribs["annotations"]

        result.update(attribs)
        return result

    def update_cds(self, line, cds):
        """Adds the fields of a gff line to an existing CDS object."""
        args = self.extract_cds_args(line)
        cds.add_indices(args['indices'])
        cds.add_phase(args['phase'])
        cds.add_identifier(args['identifier'])
        if 'score' in args:
            cds.add_score(args['score'])
        cds.sort_attributes()

    def update_exon(self, line, exon):
        """Adds the fields of a gff line to an existing Exon object."""
        args = self.extract_exon_args(line)
        exon.add_indices(args['indices'])
        exon.add_identifier(args['identifier'])
        if 'score' in args:
            exon.add_score(args['score'])
        exon.sort_attributes()

    def process_line(self, line):
        """Processes the contents of one line of a .gff file, returning True if successful.

        Args:
            line: a list of the fields
        """
        ltype = self.line_type(line)
        if ltype in ('gene', 'pseudogene', 'transposable_element_gene'):
            self.process_gene_line(line, ltype)
            return True
        elif ltype in ('mRNA', 'tRNA', 'rRNA', 'ncRNA', 'miRNA', \
                'miRNA_primary_transcript', 'snRNA', 'snoRNA', \
                'lnc_RNA', 'antisense_lncRNA', 'antisense_RNA',
                'uORF', 'pseudogenic_transcript', 'pseudogenic_tRNA', \
                'transcript_region'):
            ncrna_class = None
            if ltype in ('snRNA', 'snoRNA', 'miRNA', 'antisense_RNA', 'lnc_RNA', 'antisense_lncRNA'):
                if ltype in ('lnc_RNA', 'antisense_lncRNA'):
                    ltype = 'lncRNA'
                ltype, ncrna_class = 'ncRNA', ltype
            elif ltype == 'miRNA_primary_transcript':
                ltype = 'precursor_RNA'
            elif ltype in ('pseudogenic_transcript', 'pseudogenic_tRNA', \
                    'transcript_region'):
                ltype = 'misc_RNA'
            elif ltype == 'uORF':
                ltype = 'misc_feature'
            if ltype == 'ncRNA' and not ncrna_class:
                ncrna_class = 'other'
            self.process_rna_line(line, ltype, ncrna_class)
            return True
        elif ltype == 'CDS':
            self.process_cds_line(line)
            return True
        elif ltype in ('exon', 'pseudogenic_exon'):
            self.process_exon_line(line)
            return True
        elif ltype in ('start_codon', 'stop_codon'):
            self.process_other_feature_line(line)
            return True
        elif ltype in ('internal_eliminated_sequence', 'chromosome_breakage_sequence'):
            if ltype == 'internal_eliminated_sequence':
                ltype = 'iDNA'
            elif ltype == 'chromosome_breakage_sequence':
                ltype = 'misc_feature'
            self.process_non_gene_feature_line(line, ltype)
            return True
        else:
            self.skipped_features += 1
            return False

    def process_gene_line(self, line, gene_type):
        """Extracts arguments from a line and instantiates a Gene object."""
        kwargs = self.extract_gene_args(line)
        if not kwargs:
            return
        gene_id = kwargs['identifier']
        # noinspection PyArgumentList
        gene = Gene(**kwargs)
        if gene_type in ('pseudogene', 'transposable_element_gene'):
            gene.pseudo = True
        self.genes[gene_id] = gene

    def process_rna_line(self, line, rna_type, ncrna_class):
        """Extracts arguments from a line and instantiates an XRNA object."""
        kwargs = self.extract_mrna_args(line)
        if not kwargs:
            return
        kwargs["rna_type"] = rna_type
        kwargs["ncrna_class"] = ncrna_class
        mrna_id = kwargs['identifier']
        # noinspection PyArgumentList
        self.mrnas[mrna_id] = XRNA(**kwargs)

    def process_cds_line(self, line):
        """Extracts arguments from a line and adds them to a CDS, or makes a new one."""
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
            # noinspection PyArgumentList
            parent_mrna.cds = CDS(**kwargs)

    def process_exon_line(self, line):
        """Extracts arguments from a line and adds them to a Exon, or makes a new one."""
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
            # noinspection PyArgumentList
            parent_mrna.exon = Exon(**kwargs)

    def process_other_feature_line(self, line):
        """Extracts arguments from a line and instantiates a GenePart from them."""
        kwargs = self.extract_other_feature_args(line)
        if not kwargs:
            return
        parent_id = kwargs['parent_id']
        if parent_id not in self.mrnas:
            self.orphans.append(line)
            return
        parent_mrna = self.mrnas[parent_id]
        # noinspection PyArgumentList
        parent_mrna.other_features.append(GenePart(**kwargs))

    def process_non_gene_feature_line(self, line, feat_type):
        """Extracts arguments from a line and add them to an XRNA object without a parent Gene"""
        kwargs = self.extract_mrna_args(line)
        if not kwargs:
            return
        kwargs["rna_type"] = feat_type
        mrna_id = kwargs['identifier']
        self.mrnas[mrna_id] = XRNA(**kwargs)

    def read_file(self, reader):
        """GFFReader's public method, takes a reader, returns list of Genes,
        list of comments, list of invalid lines and list of ignored features.

        Writes comments to 'genome.comments.gff',
        invalid lines to 'genome.invalid.gff' and
        ignored features to 'genome.ignored.gff'.
        """
        # Open files to write comments, invalid and ignored entries
        comments = []
        invalid = []
        ignored = []

        # First pass, pulling out all genes and mRNAs
        #  and placing child features if possible
        for line in reader:
            if len(line) == 0 or line.startswith('#'):
                comments.append(line)
                continue
            splitlines = self.validate_line(line)
            if not splitlines:
                invalid.append(line)
            else:
                for splitline in splitlines:
                    line_added = self.process_line(splitline)
                    if not line_added:
                        ignored.append(line)

        # Second pass, placing child features which
        # preceded their parents in the first pass
        orphans = copy.deepcopy(self.orphans)
        for splitline in orphans:
            self.process_line(splitline)

        # Add mRNAs to their parent genes
        # Skip if RNA feature has a Derives_from relationshop
        for mrna in self.mrnas.values():
            if not mrna.parent_id:
                self.non_genes.append(mrna)
            else:
                parent_gene = self.genes[mrna.parent_id]
                parent_gene.mrnas.append(mrna)

        if self.skipped_features > 0:
            sys.stderr.write("Warning: skipped " + str(self.skipped_features) + " uninteresting features.\n")
        return self.genes.values(), self.non_genes, comments, invalid, ignored

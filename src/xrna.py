#!/usr/bin/env python
# coding=utf-8

import math
import re
from src.gene_part import GenePart
import src.translator as translate


def length_of_segment(index_pair):
    return math.fabs(index_pair[1] - index_pair[0]) + 1


class XRNA(object):
    def __init__(self, identifier, indices, parent_id, source='', seq_name='', name='', strand='+',
                 annotations=None, rna_type="mRNA", ncrna_class=None, derives_from=None):
        self.rna_type = rna_type
        self.ncrna_class = ncrna_class
        self.derives_from = derives_from
        self.identifier = identifier
        self.indices = indices
        self.parent_id = parent_id
        self.strand = strand
        self.exon = None
        self.cds = None
        self.other_features = []
        self.annotations = {} if annotations is None else annotations
        self.death_flagged = False
        self.source = source
        self.seq_name = seq_name
        self.name = name

    def __str__(self):
        """Returns string representation of the RNA.

        String contains the RNA's identifier and the number of features it contains.
        """
        result = self.rna_type + " (ID=" + str(self.identifier) + ") containing "
        if self.exon:
            result += "Exon, "
        if self.cds:
            result += "CDS "
        if len(self.other_features) > 0:
            result += "and " + str(len(self.other_features))
            result += " other features"
        return result

    def add_annotation(self, key, value, feat_type=None):
        """Adds an annotation to the RNA.

        Args:
            key: the type of annotation
            value: the annotation itself
            feat_type: type of RNA or sub feature
        """
        if not feat_type:
            if key not in self.annotations:
                self.annotations[key] = []

            self.annotations[key].append(value)
        else:
            if feat_type not in self.annotations:
                self.annotations[feat_type] = dict()

            if key not in self.annotations[feat_type]:
                self.annotations[feat_type][key] = []

            self.annotations[feat_type][key].append(value)

    def length(self):
        """Returns the length of the RNA."""
        return length_of_segment(self.indices)

    def adjust_indices(self, n, start_index=1):
        """Increments indices of RNA and its child features.

        Optionally, only indices occurring after start_index are incremented.

        Args:
            n: integer by which to increment indices
            start_index: optional coordinate before which no indices will be changed.
        """
        if self.indices[0] > start_index:
            self.indices = [i + n for i in self.indices]
        elif self.indices[1] > start_index:
            self.indices[1] += n
        if self.exon:
            self.exon.adjust_indices(n, start_index)
        if self.cds:
            self.cds.adjust_indices(n, start_index)
        for feature in self.other_features:
            feature.adjust_indices(n, start_index)

    def number_of_gagflags(self):
        """Returns the number of flagged features contained by RNA.

        Multiple flags on a given CDS or Exon are ignored, so the
        possible return values are 0, 1 or 2.
        """
        total = 0
        if self.cds and self.cds.gagflagged():
            total += 1
        if self.exon and self.exon.gagflagged():
            total += 1
        return total

    def create_start_and_stop_if_necessary(self, seq_object, strand):
        """Inspects child CDS and creates start/stop codons if appropriate.

        This is accomplished by examining the first and last three nucleotides
        of the CDS and comparing them to start/stop codon sequences. This will
        override any start_codon or stop_codon entries in the original GFF file.

        Args:
            seq_object: the actual sequence containing the RNA
            strand: either '+' or '-'
        """
        # TODO I'd rather pass seq.bases than the object itself, since
        # the object owns this mrna...
        if not self.cds:
            return
        seq = self.cds.extract_sequence(seq_object, strand)
        if translate.has_start_codon(seq):
            indices = self.cds.get_start_indices(strand)
            self.add_start_codon(indices)
        if translate.has_stop_codon(seq):
            indices = self.cds.get_stop_indices(strand)
            self.add_stop_codon(indices)

    def indices_intersect_mrna(self, indices):
        """Returns a boolan indicating whether a pair of indices overlap the RNA"""
        if len(indices) != 2:
            return False
        begin = indices[0]
        end = indices[1]
        self_start = self.indices[0]
        self_end = self.indices[1]
        # mrna contains beginning of indices
        if self_start <= begin <= self_end:
            return True
        # mrna contains end of indices
        elif self_start <= end <= self_end:
            return True
        # indices contain entire mrna
        elif begin <= self_start and end >= self_end:
            return True
        else:
            return False

    def add_other_feature(self, feature):
        """Adds a feature to MRNA.other_features list"""
        self.other_features.append(feature)

    def add_start_codon(self, indices):
        """Adds a start_codon GenePart to MRNA.other_features"""
        # TODO figure out naming scheme...
        start_id = self.identifier + ":start"
        start_parent_id = self.identifier
        start = GenePart(feature_type='start_codon', identifier=start_id,
                         indices=indices, parent_id=start_parent_id, strand=self.strand)
        self.add_other_feature(start)

    def add_stop_codon(self, indices):
        """Adds a stop_codon GenePart to MRNA.other_features"""
        stop_id = self.identifier + ":stop"
        stop_parent_id = self.identifier
        stop = GenePart(feature_type='stop_codon', identifier=stop_id,
                        indices=indices, parent_id=stop_parent_id, strand=self.strand)
        self.add_other_feature(stop)

    def has_start(self):
        """Returns a boolean indicating whether the RNA contains a 'start_codon' feature"""
        if self.rna_type != 'mRNA': return True
        for feature in self.other_features:
            if feature.feature_type == 'start_codon':
                return True
        return False

    def has_stop(self):
        """Returns a boolean indicating whether the RNA contains a 'stop_codon' feature"""
        if self.rna_type != 'mRNA': return True
        for feature in self.other_features:
            if feature.feature_type == 'stop_codon':
                return True
        return False

    def indices_intersect_cds(self, indices):
        """Returns a boolean indicating whether given indices overlap RNA's CDS"""
        if not self.cds:
            return False
        else:
            return self.cds.indices_intersect_cds(indices)

    def cds_to_gff(self, seq_id, source):
        """Returns a string containing RNA's child CDS in .gff format."""
        if self.cds:
            return self.cds.to_gff(seq_id, source)
        else:
            return ""

    def cds_to_tbl(self):
        """Returns a string containing RNA's child CDS in .tbl format."""
        if self.cds:
            has_start = self.has_start()
            has_stop = self.has_stop()
            return self.cds.to_tbl(has_start, has_stop)
        else:
            return ""

    def to_gff(self):
        """Returns a string of RNA and child features in .gff format."""
        result = self.seq_name + "\t" + self.source + "\t" + self.rna_type + "\t"
        result += str(self.indices[0]) + "\t" + str(self.indices[1]) + "\t"
        result += "." + "\t" + self.strand + "\t" + "." + "\t"
        result += "ID=" + str(self.identifier)
        result += ";Parent=" + str(self.parent_id)
        xrna_annotations = self.annotations[self.rna_type] \
                if self.rna_type in self.annotations.keys() \
                else self.annotations
        for key in xrna_annotations.keys():
            result += ';' + key + "="
            result += ','.join(xrna_annotations[key])
        result += '\n'
        if self.exon:
            result += self.exon.to_gff(self.seq_name, self.source)
        if self.cds:
            result += self.cds.to_gff(self.seq_name, self.source)
        for other in self.other_features:
            result += other.to_gff(self.seq_name, self.source)
        return result

    def to_tbl(self, gc_tag='ncbi', txid_format='suffix'):
        """Returns a string of RNA and child features in .tbl format."""
        has_start = self.has_start()
        has_stop = self.has_stop()
        output = ""

        # set the transcript_id based on specified format
        transcript_id = "{0}_mrna".format(self.identifier) \
            if txid_format == 'suffix' \
            else 'mRNA.{0}'.format(self.identifier)

        if self.exon or self.rna_type in ('precursor_RNA', 'misc_feature', 'iDNA') \
                or self.ncrna_class == 'miRNA':
            if self.exon:
                output += self.exon.to_tbl(has_start, has_stop, self.rna_type)
            else:
                (start, stop) = (self.indices[0], self.indices[1]) if self.strand == '+' else \
                    (self.indices[1], self.indices[0])
                output += str(start)+'\t'+str(stop)+'\t'+self.rna_type+'\n'

            # Write the annotations
            if not self.ncrna_class and not self.annotations_contain_product(feat_type=self.rna_type) \
                    and not self.rna_type in ('misc_RNA', 'misc_feature', 'iDNA'):
                output += self.tbl_line("product", "hypothetical protein")

            if self.rna_type == "mRNA":
                output += self.tbl_line("protein_id", "gnl|{0}|{1}".format(gc_tag, self.identifier))

            if self.rna_type in ('misc_RNA', 'misc_feature', 'iDNA'):
                if self.rna_type == 'iDNA':
                    identifier = "{0} ({1})".format("Internal Eliminated Sequence", self.identifier)
                elif self.rna_type == 'misc_feature' and re.search("[0-9]+[LR]-[0-9]+", self.identifier):
                    identifier = "{0} (Cbs {1})".format("Chromosome breakage sequence", self.identifier)
                else:
                    identifier = self.identifier
                output += self.tbl_line("note", identifier)
            elif not self.ncrna_class == 'miRNA':
                output += self.tbl_line("transcript_id", "gnl|{0}|{1}".format(gc_tag, transcript_id))

            if self.ncrna_class:
                output += self.tbl_line("ncRNA_class", self.ncrna_class)
                if not self.annotations_contain_product(feat_type=self.rna_type):
                    output += self.tbl_line("product", "other RNA")

            xrna_annotations = self.annotations[self.rna_type] \
                    if self.rna_type in self.annotations.keys() \
                    else self.annotations

            for key in xrna_annotations.keys():
                for value in xrna_annotations[key]:
                    output += self.tbl_line(key, value)
        if self.cds:
            output += self.cds.to_tbl(has_start, has_stop)
            # Write the annotations
            cds_annotations = self.annotations['CDS'] \
                    if 'CDS' in self.annotations.keys() \
                    else self.annotations

            for key in sorted(cds_annotations.keys()):
                if key == 'protein_id':
                    continue
                for value in cds_annotations[key]:
                    output += self.tbl_line(key, value)

            if not self.annotations_contain_product(feat_type="CDS"):
                output += self.tbl_line("product", "hypothetical protein")

            if self.annotations_contain_protein_id(feat_type="CDS"):
                output += self.tbl_line("protein_id", \
                    "gnl|{0}|{1}|gb|{2}".format(gc_tag, self.identifier, cds_annotations['protein_id'][0]))
            else:
                output += self.tbl_line("protein_id", "gnl|{0}|{1}".format(gc_tag, self.identifier))
            output += self.tbl_line("transcript_id", "gnl|{0}|{1}".format(gc_tag, transcript_id))
        return output

    def tbl_line(self, *args, **kwargs):
        return "\t\t\t{0}\t{1}\n".format(args[0], args[1]) \
            if len(args) == 2 else \
            "\t\t\t{0}\n".format(args[0])

    ## STATS STUFF ##

    def get_longest_exon(self):
        """Returns length of longest exon contained on RNA."""
        if not self.exon:
            return 0
        longest = 0
        for index_pair in self.exon.indices:
            if length_of_segment(index_pair) > longest:
                longest = length_of_segment(index_pair)
        return longest

    def get_shortest_exon(self):
        """Returns length of shortest exon contained on RNA."""
        if not self.exon:
            return 0
        shortest = None
        for index_pair in self.exon.indices:
            length = length_of_segment(index_pair)
            if length == 0:
                continue
            if shortest is None or length_of_segment(index_pair) < shortest:
                shortest = length
        if shortest is None:
            return 0
        return shortest

    def get_total_exon_length(self):
        """Returns sum of all child exon lengths."""
        if not self.exon:
            return 0

        total = 0
        for index_pair in self.exon.indices:
            total += length_of_segment(index_pair)
        return total

    def get_num_exons(self):
        """Returns number of exons contained on RNA."""
        if self.exon:
            return len(self.exon.indices)
        else:
            return 0

    def get_longest_intron(self):
        """Returns length of longest intron contained on RNA."""
        if not self.exon:
            return 0
        longest = 0
        last_end = 0
        for index_pair in self.exon.indices:
            if last_end != 0:
                this_intron = index_pair[0] - last_end - 1
                if this_intron > longest:
                    longest = this_intron
            last_end = index_pair[1]
        return longest

    def get_shortest_intron(self):
        """Returns length of shortest intron contained on RNA."""
        if not self.exon:
            return 0
        shortest = None
        last_end = 0
        for index_pair in self.exon.indices:
            if last_end != 0:
                this_intron = index_pair[0] - last_end - 1
                if this_intron == 0:
                    continue
                if this_intron < 0:
                    raise Exception("Intron with negative length on " + self.identifier)
                if shortest is None or this_intron < shortest:
                    shortest = this_intron
            last_end = index_pair[1]
        if shortest is None:
            return 0
        return shortest

    def get_total_intron_length(self):
        """Returns sum of lengths of all introns contained on RNA."""
        if not self.exon:
            return 0
        total = 0
        last_end = 0
        for index_pair in self.exon.indices:
            if last_end != 0:
                total += abs(index_pair[0] - last_end) + 1
            last_end = index_pair[1]
        return total

    def get_num_introns(self):
        """Returns number of introns contained on RNA."""
        if self.exon:
            return len(self.exon.indices) - 1
        else:
            return 0

    def annotations_contain_product(self, feat_type=None):
        if feat_type not in self.annotations.keys():
            return 'product' in self.annotations.keys()

        return 'product' in self.annotations[feat_type].keys()

    def annotations_contain_protein_id(self, feat_type=None):
        if feat_type not in self.annotations.keys():
            return 'protein_id' in self.annotations.keys()

        return 'protein_id' in self.annotations[feat_type].keys()

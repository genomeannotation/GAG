#!/usr/bin/env python

import math
import sys
from feature_tbl_entry import FeatureTblEntry
from src.gene_part import GenePart
from src.translate import *

def length_of_segment(index_pair):
    return math.fabs(index_pair[1] - index_pair[0]) + 1

class MRNA:

    def __init__(self, identifier, name, indices, parent_id):
        self.identifier = identifier
        self.name = name
        self.indices = indices
        self.parent_id = parent_id
        self.exon = None
        self.cds = None
        self.other_features = []

    def __str__(self):
        result = "mRNA (ID=" + str(self.identifier) + ", Name="
        result += self.name + ") containing "
        if self.exon:
            result += "Exon, "
        if self.cds:
            result += "CDS "
        if len(self.other_features) > 0:
            result += "and " + str(len(self.other_features)) 
            result += " other features"
        return result

    def length(self):
        return length_of_segment(self.indices)

    def is_maker_mrna(self):
        return 'maker' in self.name

    def adjust_indices(self, n):
        self.indices = [i + n for i in self.indices]
        if self.exon:
            self.exon.adjust_indices(n)
        if self.cds:
            self.cds.adjust_indices(n)
        for feature in self.other_features:
            feature.adjust_indices(n)

    def adjust_phase(self):
        if self.cds:
            self.cds.adjust_phase()

    def clean_up_indices(self):
        if self.indices[1] < 1:
            self.indices[0] = 0
            self.indices[1] = 0
        elif self.indices[0] < 1:
            self.indices[0] = 1
        if self.cds:
            self.cds.clean_up_indices()
        if self.exon:
            self.exon.clean_up_indices()
        for feat in self.other_features:
            feat.clean_up_indices()

    def create_start_and_stop_if_necessary(self, fasta, seq_name, phase):
        if not self.cds:
            return
        seq = self.cds.extract_sequence(fasta, seq_name, phase)
        if has_start_codon(seq):
            indices = self.cds.get_start_indices(phase)
            self.add_start_codon(indices)
        if has_stop_codon(seq):
            indices = self.cds.get_stop_indices(phase)
            self.add_stop_codon(indices)


    def remove_first_cds_segment_if_shorter_than(self, min_length):
        if self.cds:
            if length_of_segment(self.cds.indices[0]) < min_length:
                self.cds.indices = self.cds.indices[1:]

    def trim_end(self, endindex):
        if self.indices[0] > endindex:
            self.indices[0] = 0
            self.indices[1] = 0
        elif self.indices[1] > endindex:
            self.indices[1] = endindex
            if self.cds:
                self.cds.trim_end(endindex)
            if self.exon:
                self.exon.trim_end(endindex)
            for feat in self.other_features:
                feat.trim_end(endindex)

    def remove_invalid_features(self):
        if self.cds:
            self.cds.remove_trimmed_segments()
            ind = self.cds.indices
            if len(ind) == 0:
                self.cds = None
        if self.exon:
            self.exon.remove_trimmed_segments()
            if len(self.exon.indices) == 0:
                self.exon = None
        invalid_features = []
        for i, feat in enumerate(self.other_features):
            self.other_features[i].remove_trimmed_segments()
            if len(feat.indices) == 0:
                invalid_features.append(i)
            elif feat.feature_type == 'start_codon' or feat.feature_type == 'stop_codon':
                if not feat.valid_codon():
                    invalid_features.append(i)
        for j in reversed(invalid_features):
            self.other_features.pop(j)


    def set_exon(self, exon):
        self.exon = exon

    def set_cds(self, cds):
        self.cds = cds

    def add_other_feature(self, feature):
        self.other_features.append(feature)

    def add_start_codon(self, begin_index):
        start_id = self.identifier + 1000000
        start_name = self.name + ':start'
        start_indices = [begin_index, begin_index+2]
        start_parent_id = self.identifier
        start = GenePart(feature_type='start_codon', identifier=start_id, name=start_name, indices=start_indices, parent_id=start_parent_id)
        self.add_other_feature(start)

    def add_stop_codon(self, end_index):
        stop_id = self.identifier + 1000001
        stop_name = self.name + ':stop'
        stop_indices = [end_index-2, end_index]
        stop_parent_id = self.identifier
        stop = GenePart(feature_type='stop_codon', identifier=stop_id, name=stop_name, indices=stop_indices, parent_id=stop_parent_id)
        self.add_other_feature(stop)

    def get_cds_indices(self):
        if self.cds:
            result = []
            for i, index_pair in enumerate(self.cds.indices):
                pair = self.cds.indices[i]
                new_start = pair[0] + self.cds.get_phase(i)
                result.append([new_start, pair[1]])
            return result
        else:
            return None

    def length_of_shortest_cds_segment(self):
        return self.cds.length_of_shortest_segment()

    def has_start(self):
        for feature in self.other_features:
            if feature.feature_type is 'start_codon':
                return True
        return False

    def has_stop(self):
        for feature in self.other_features:
            if feature.feature_type is 'stop_codon':
                return True
        return False

    def to_gff(self, seq_name, source, strand):
        result = seq_name + "\t" + source + "\t" + "mRNA" + "\t"
        result += str(self.indices[0]) + "\t" + str(self.indices[1]) + "\t"
        result += "." + "\t" + strand + "\t" + "." + "\t"
        result += "ID=" + str(self.identifier) + ";Name=" + self.name
        result += ";Parent=" + str(self.parent_id) + "\n"
        if self.exon:
            result += self.exon.to_gff(seq_name, source, strand)
        if self.cds:
            result += self.cds.to_gff(seq_name, source, strand)
        for other in self.other_features:
            result += other.to_gff(seq_name, source, strand)
        return result

    def to_tbl_entries(self, strand):
        entries = []

        has_start = False
        has_stop = False

        for other in self.other_features:
            if other.feature_type == 'start_codon':
                has_start = True
            elif other.feature_type == 'stop_codon':
                has_stop = True

        phase = 0
        if self.cds != None:
            phase = self.cds.get_phase(0)
            cdsEntry = FeatureTblEntry()
            cdsEntry.set_type("CDS")
            cdsEntry.set_name(self.name)
            for coord in self.cds.indices:
                cdsEntry.add_coordinates(coord[0], coord[1])
            cdsEntry.set_strand(strand)
            cdsEntry.set_phase(phase)
            cdsEntry.set_partial_info(has_start, has_stop)
            entries.append(cdsEntry)

        if self.exon != None:
            exonEntry = FeatureTblEntry()
            exonEntry.set_type("mRNA")
            exonEntry.set_name(self.name)
            for coord in self.exon.indices:
                exonEntry.add_coordinates(coord[0], coord[1])
            exonEntry.set_strand(strand)
            exonEntry.set_phase(phase)
            exonEntry.set_partial_info(has_start, has_stop)
            entries.append(exonEntry)

        return entries

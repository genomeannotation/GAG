#!/usr/bin/env python

import math
import sys

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
        if self.exon or self.cds:
            result += "and "
        if len(self.other_features) > 0:
            result += str(len(self.other_features)) + " other features"
        return result

    def length(self):
        return length_of_segment(self.indices)

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
            elif feat.feature_type == 'start_codon' or 'stop_codon':
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

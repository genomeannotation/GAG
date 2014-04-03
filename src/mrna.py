#!/usr/bin/env python

import math
from src.feature_tbl_entry import FeatureTblEntry
from src.gene_part import GenePart
import src.translate as translate

def length_of_segment(index_pair):
    return math.fabs(index_pair[1] - index_pair[0]) + 1

class MRNA:

    def __init__(self, identifier, indices, parent_id, annotations=[]):
        self.identifier = identifier
        self.indices = indices
        self.parent_id = parent_id
        self.exon = None
        self.cds = None
        self.other_features = []
        self.annotations = annotations
        self.death_flagged = False

    def __str__(self):
        result = "mRNA (ID=" + str(self.identifier) + ") containing "
        if self.exon:
            result += "Exon, "
        if self.cds:
            result += "CDS "
        if len(self.other_features) > 0:
            result += "and " + str(len(self.other_features)) 
            result += " other features"
        return result
        
    def add_annotation(self, key, value):
        self.annotations.append([key, value])

    def length(self):
        return length_of_segment(self.indices)

    def is_maker_mrna(self):
        return 'maker' in self.identifier

    def adjust_indices(self, n, start_index=1):
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

    def create_start_and_stop_if_necessary(self, seq_object, phase):
        # TODO I'd rather pass seq.bases than the object itself, since
        # the object owns this mrna...
        if not self.cds:
            return
        seq = self.cds.extract_sequence(seq_object, phase)
        if translate.has_start_codon(seq):
            indices = self.cds.get_start_indices(phase)
            self.add_start_codon(indices)
        if translate.has_stop_codon(seq):
            indices = self.cds.get_stop_indices(phase)
            self.add_stop_codon(indices)


    def remove_first_cds_segment_if_shorter_than(self, min_length):
        if self.cds:
            if length_of_segment(self.cds.indices[0]) < min_length:
                self.cds.indices = self.cds.indices[1:]
                self.cds.identifier = self.cds.identifier[1:]
                self.cds.phase = self.cds.phase[1:]

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
            elif feat.feature_type == 'start_codon' \
                    or feat.feature_type == 'stop_codon':
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

    def add_start_codon(self, indices):
        # TODO figure out naming scheme...
        start_id = self.identifier + ":start"
        start_parent_id = self.identifier
        start = GenePart(feature_type='start_codon', identifier=start_id, \
                indices=indices, parent_id=start_parent_id)
        self.add_other_feature(start)

    def add_stop_codon(self, indices):
        stop_id = self.identifier + ":stop"
        stop_parent_id = self.identifier
        stop = GenePart(feature_type='stop_codon', identifier=stop_id, \
                indices=indices, parent_id=stop_parent_id)
        self.add_other_feature(stop)

    def length_of_shortest_cds_segment(self):
        return self.cds.length_of_shortest_segment()

    def has_start(self):
        for feature in self.other_features:
            if feature.feature_type == 'start_codon':
                return True
        return False

    def has_stop(self):
        for feature in self.other_features:
            if feature.feature_type == 'stop_codon':
                return True
        return False

    def get_longest_exon(self):
        longest = 0
        for index_pair in self.exon.indices:
            if length_of_segment(index_pair) > longest:
                longest = length_of_segment(index_pair)
        return longest

    def get_shortest_exon(self):
        shortest = 0
        for index_pair in self.exon.indices:
            length = length_of_segment(index_pair)
            if shortest == 0 or length_of_segment(index_pair) < shortest:
                shortest = length_of_segment(index_pair)
        return shortest

    def to_gff(self, seq_name, source, strand, death_flagged_stuff=False):
        if not death_flagged_stuff and self.death_flagged:
            return ""
    
        result = seq_name + "\t" + source + "\t" + "mRNA" + "\t"
        result += str(self.indices[0]) + "\t" + str(self.indices[1]) + "\t"
        result += "." + "\t" + strand + "\t" + "." + "\t"
        result += "ID=" + str(self.identifier)
        result += ";Parent=" + str(self.parent_id)
        for annot in self.annotations:
            result += ';'+annot[0]+'='+annot[1]
        result += '\n'
        if self.exon:
            result += self.exon.to_gff(seq_name, source, strand)
        if self.cds:
            result += self.cds.to_gff(seq_name, source, strand)
        for other in self.other_features:
            result += other.to_gff(seq_name, source, strand)
        return result

    def to_tbl_entries(self, annotator, strand):
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
            cdsEntry.set_name(self.identifier)
            for coord in self.cds.indices:
                cdsEntry.add_coordinates(coord[0], coord[1])
            cdsEntry.set_strand(strand)
            cdsEntry.set_phase(phase)
            cdsEntry.set_partial_info(has_start, has_stop)
            annotator.annotate_cds(cdsEntry)
            entries.append(cdsEntry)

        if self.exon != None:
            exonEntry = FeatureTblEntry()
            exonEntry.set_type("mRNA")
            exonEntry.set_name(self.identifier)
            for coord in self.exon.indices:
                exonEntry.add_coordinates(coord[0], coord[1])
            exonEntry.set_strand(strand)
            exonEntry.set_phase(phase)
            exonEntry.set_partial_info(has_start, has_stop)
            annotator.annotate_mrna(exonEntry)
            entries.append(exonEntry)

        return entries

    def to_tbl(self, strand):
        if self.death_flagged:
            return ""
    
        has_start = self.has_start()
        has_stop = self.has_stop()
        output = ""
        if self.exon:
            output += self.exon.to_tbl(strand, has_start, has_stop)
        if self.cds:
            output += self.cds.to_tbl(strand, has_start, has_stop)
        return output

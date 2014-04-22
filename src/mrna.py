#!/usr/bin/env python

import math
from src.gene_part import GenePart
import src.translator as translate

def length_of_segment(index_pair):
    return math.fabs(index_pair[1] - index_pair[0]) + 1

class MRNA:

    def __init__(self, identifier, indices, parent_id, strand='+', annotations=None):
        self.identifier = identifier
        self.indices = indices
        self.parent_id = parent_id
        self.strand = strand
        self.exon = None
        self.cds = None
        self.other_features = []
        if not annotations:
            self.annotations = []
        else:
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

    def indices_intersect_mrna(self, indices):
        if len(indices) != 2:
            return False
        begin = indices[0]
        end = indices[1]
        self_start = self.indices[0]
        self_end = self.indices[1]
        # mrna contains beginning of indices
        if self_start <= begin and self_end >= begin:
            return True
        # mrna contains end of indices
        elif self_start <= end and self_end >= end:
            return True
        # indices contain entire mrna
        elif begin <= self_start and end >= self_end:
            return True
        else:
            return False

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

    def indices_intersect_cds(self, indices):
        if not self.cds:
            return False
        else:
            return self.cds.indices_intersect_cds(indices)

    def cds_to_gff(self, seq_id, source):
        if self.cds:
            return self.cds.to_gff(seq_id, source)
        else:
            return ""

    def to_gff(self, seq_name, source):
        result = seq_name + "\t" + source + "\t" + "mRNA" + "\t"
        result += str(self.indices[0]) + "\t" + str(self.indices[1]) + "\t"
        result += "." + "\t" + self.strand + "\t" + "." + "\t"
        result += "ID=" + str(self.identifier)
        result += ";Parent=" + str(self.parent_id)
        for annot in self.annotations:
            result += ';'+annot[0]+'='+annot[1]
        result += '\n'
        if self.exon:
            result += self.exon.to_gff(seq_name, source)
        if self.cds:
            result += self.cds.to_gff(seq_name, source)
        for other in self.other_features:
            result += other.to_gff(seq_name, source)
        return result

    def to_tbl(self, strand):
        has_start = self.has_start()
        has_stop = self.has_stop()
        output = ""
        if self.exon:
            output += self.exon.to_tbl(strand, has_start, has_stop)
        if self.cds:
            output += self.cds.to_tbl(strand, has_start, has_stop)
        return output

    ## STATS STUFF ##

    def get_longest_exon(self):
        if not self.exon:
            return 0
    
        longest = 0
        for index_pair in self.exon.indices:
            if length_of_segment(index_pair) > longest:
                longest = length_of_segment(index_pair)
        return longest

    def get_shortest_exon(self):
        if not self.exon:
            return 0
        shortest = 0
        for index_pair in self.exon.indices:
            length = length_of_segment(index_pair)
            if shortest == 0 or length_of_segment(index_pair) < shortest:
                shortest = length
        return shortest

    def get_total_exon_length(self):
        if not self.exon:
            return 0
    
        total = 0
        for index_pair in self.exon.indices:
            total += length_of_segment(index_pair)
        return total

    def get_num_exons(self):
        if self.exon:
            return len(self.exon.indices)
        else:
            return 0

    def get_longest_intron(self):
        if not self.exon:
            return 0
        longest = 0
        last_end = 0
        for index_pair in self.exon.indices:
            if last_end != 0:
                this_intron = abs(index_pair[0] - last_end) + 1
                if this_intron > longest:
                    longest = this_intron
            last_end = index_pair[1]
        return longest

    def get_shortest_intron(self):
        if not self.exon:
            return 0
        shortest = 0
        last_end = 0
        for index_pair in self.exon.indices:
            if last_end != 0:
                this_intron = abs(index_pair[0] - last_end) + 1
                if shortest == 0 or this_intron < shortest:
                    shortest = this_intron
            last_end = index_pair[1]
        return shortest

    def get_total_intron_length(self):
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
        if self.exon:
            return len(self.exon.indices) - 1
        else:
            return 0

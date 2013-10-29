#!/usr/bin/env python

import math

def length_of_segment(index_pair):
    return math.fabs(index_pair[1] - index_pair[0]) + 1

def adjust_index_pair(index_pair, n):
    return [i + n for i in index_pair]

class GenePart:
    def __init__(self, feature_type=None, id=[], name=[], indices=[], score=[], parent_id=None):
        self.feature_type = feature_type
        self.id = id
        self.name = name
        self.indices = indices
        self.score = score
        self.parent_id = parent_id

    def length(self):
        length = 0
        for index_pair in self.indices:
            length += length_of_segment(index_pair)
        return length

    # used by .to_gff
    def get_score(self, i):
        if self.score and len(self.score) > i:
            return self.score[i]
        else:
            return "."

    # used by .to_gff
    # (CDS overrides this method)
    def get_phase(self, i):
        return "."

    def adjust_indices(self, n):
        self.indices = [adjust_index_pair(pair, n) for pair in self.indices]

    def length_of_shortest_segment(self):
        if len(self.indices) == 0:
            return None
        min_length = length_of_segment(self.indices[0])
        if len(self.indices) == 1:
            return min_length
        else:
            for index_pair in self.indices:
                if length_of_segment(index_pair) < min_length:
                    min_length = length_of_segment(index_pair)
        return min_length

    def generate_attribute_entry(self, i):
        if len(self.id) <= i or self.parent_id is None:
            return None
        entry = "ID=" + str(self.id[i]) + ";"
        if len(self.name) > i:
            entry += "Name=" + str(self.name[i]) + ";"
        entry += "Parent=" + str(self.parent_id) + "\n"
        return entry

    def to_gff(self, seq_name, source, strand):
        result = ""
        for i in range(len(self.indices)):
            result += seq_name + "\t" + source + "\t"
            result += self.feature_type + "\t" + str(self.indices[i][0])
            result += "\t" + str(self.indices[i][1]) + "\t"
            result += str(self.get_score(i)) + "\t" + strand + "\t"
            result += str(self.get_phase(i)) + "\t"
            result += self.generate_attribute_entry(i)
        return result


class CDS(GenePart):

    def __init__(self, id=[], name=[], indices=[], score=[], phase=[], parent_id=None):
        GenePart.__init__(self, feature_type='CDS', id=id, name=name, indices=indices, score=score, parent_id=parent_id)
        self.phase = phase 
        self.annotations = []

    def get_phase(self, i):
        if self.phase and len(self.phase) > i:
            return self.phase[i]
        else:
            return "."


class Exon(GenePart):

    def __init__(self, id=[], name=[], indices=[], score=[], parent_id=None):
        GenePart.__init__(self, feature_type='exon', id=id, name=name, indices=indices, score=score, parent_id=parent_id)
        self.annotations = []


class MRNA:

    def __init__(self, id, name, indices, parent_id):
        self.id = id
        self.name = name
        self.indices = indices
        self.parent_id = parent_id
        self.exon = None
        self.cds = None
        self.other_features = []

    def length(self):
        return length_of_segment(self.indices)

    def adjust_indices(self, n):
        self.indices = [i + n for i in self.indices]

    def set_exon(self, exon):
        self.exon = exon

    def set_cds(self, cds):
        self.cds = cds

    def add_other_feature(self, feature):
        self.other_features.append(feature)

    def length_of_shortest_cds_segment(self):
        return self.cds.length_of_shortest_segment()

    def to_gff(self, seq_name, source, strand):
        result = seq_name + "\t" + source + "\t" + "mRNA" + "\t"
        result += str(self.indices[0]) + "\t" + str(self.indices[1]) + "\t"
        result += "." + "\t" + strand + "\t" + "." + "\t"
        result += "ID=" + str(self.id) + ";Name=" + self.name
        result += ";Parent=" + str(self.parent_id) + "\n"
        result += self.exon.to_gff(seq_name, source, strand)
        result += self.cds.to_gff(seq_name, source, strand)
        for other in self.other_features:
            result += other.to_gff(seq_name, source, strand)
        return result


class Gene:

    def __init__(self, seq_name, source, indices, strand, id, name):
        self.seq_name = seq_name
        self.source = source
        self.indices = indices
        self.strand = strand
        self.id = id
        self.name = name
        self.mrnas = []

    def length(self):
        return length_of_segment(self.indices)

    def add_mrna(self, mrna):
        self.mrnas.append(mrna)

    def length_of_shortest_cds_segment(self):
        min_length = self.mrnas[0].length_of_shortest_cds_segment()
        if len(self.mrnas) == 1:
            return min_length
        else:
            for mrna in self.mrnas:
                if mrna.length_of_shortest_cds_segment() < min_length:
                    min_length = mrna.length_of_shortest_cds_segment()
        return min_length

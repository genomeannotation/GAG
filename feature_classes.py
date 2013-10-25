#!/usr/bin/env python

import math

def length_of_segment(index_pair):
    return math.fabs(index_pair[1] - index_pair[0]) + 1

def adjust_index_pair(index_pair, n):
    return [i + n for i in index_pair]

class CDS:


    def __init__(self, ids, names, indices, frames, parent_id):
        self.ids = ids
        self.names = names
        self.indices = indices
        self.frames = frames
        self.parent_id = parent_id
        self.annotations = []

    def length_of_shortest_segment(self):
        min_length = length_of_segment(self.indices[0])
        if len(self.indices) == 1:
            return min_length
        else:
            for index_pair in self.indices:
                if length_of_segment(index_pair) < min_length:
                    min_length = length_of_segment(index_pair)
        return min_length

    def length(self):
        length = 0
        for index_pair in self.indices:
            length += length_of_segment(index_pair)
        return length

    def adjust_indices(self, n):
        self.indices = [adjust_index_pair(pair, n) for pair in self.indices]

    def to_gff(self, seq_name, source, strand):
        result = ""
        for i in range(len(self.indices)):
            result += seq_name + "\t" + source + "\t" + "CDS" + "\t"
            result += str(self.indices[i][0]) + "\t" + str(self.indices[i][1]) + "\t"
            result += "." + "\t" + strand + "\t" + str(self.frames[i]) + "\t"
            result += "ID=" + str(self.ids[i]) + ";Name=" + self.names[i]
            result += ";Parent=" + str(self.parent_id) + "\n"
        return result


class Exon:

    def __init__(self, ids, names, indices, scores, parent_id):
        self.ids = ids
        self.names = names
        self.indices = indices
        self.scores = scores
        self.parent_id = parent_id
        self.annotations = []

    def length(self):
        length = 0
        for index_pair in self.indices:
            length += length_of_segment(index_pair)
        return length

    def adjust_indices(self, n):
        self.indices = [adjust_index_pair(pair, n) for pair in self.indices]

    def to_gff(self, seq_name, source, strand):
        result = ""
        for i in range(len(self.indices)):
            result += seq_name + "\t" + source + "\t" + "exon" + "\t"
            result += str(self.indices[i][0]) + "\t" + str(self.indices[i][1]) + "\t"
            result += str(self.scores[i]) + "\t" + strand + "\t" + "." + "\t"
            result += "ID=" + str(self.ids[i]) + ";Name=" + self.names[i]
            result += ";Parent=" + str(self.parent_id) + "\n"
        return result


class OtherFeature:

    def __init__(self, feature_type, indices, id, name, parent_id):
        self.feature_type = feature_type
        self.indices = indices
        self.id = id
        self.name = name
        self.parent_id = parent_id

    def length(self):
        return length_of_segment(self.indices) 

    def adjust_indices(self, n):
        self.indices = [i + n for i in self.indices]

    def to_gff(self, seq_name, source, strand):
        result = seq_name + "\t" + source + "\t" + self.feature_type + "\t"
        result += str(self.indices[0]) + "\t" + str(self.indices[1]) + "\t"
        result += "." + "\t" + strand + "\t" + "." 
        result += "ID=" + str(self.id) + ";Name=" + self.name + ";"
        result += "Parent=" + str(self.parent_id) + "\n"
        return result


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





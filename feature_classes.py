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




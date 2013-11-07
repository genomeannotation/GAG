#!/usr/bin/env python

import math
import sys

def length_of_segment(index_pair):
    return math.fabs(index_pair[1] - index_pair[0]) + 1

class Gene:

    def __init__(self, seq_name, source, indices, strand, identifier, name, score=None):
        self.seq_name = seq_name
        self.source = source
        self.indices = indices
        self.score = score
        self.strand = strand
        self.identifier = identifier
        self.name = name
        self.mrnas = []

    def __str__(self):
        result = "Gene (ID=" + str(self.identifier) + ", Name="
        result += self.name + ", seq_name=" + self.seq_name
        result += ") containing " + str(len(self.mrnas))
        result += " mrnas"
        return result

    def length(self):
        return length_of_segment(self.indices)

    def get_score(self):
        if self.score:
            return self.score
        else:
            return '.'

    def add_mrna(self, mrna):
        self.mrnas.append(mrna)

    def trim_end(self, endindex):
        if self.indices[0] > endindex:
            self.indices[0] = 0
            self.indices[1] = 0
        elif self.indices[1] > endindex:
            self.indices[1] = endindex
            for mrna in self.mrnas:
                mrna.trim_end(endindex)

    # beginindex is the new start index of sequence
    def trim_begin(self, beginindex):
        self.adjust_indices(-beginindex + 1)

    def clean_up_indices(self):
        if self.indices[1] < 1:
            self.indices[0] = 0
            self.indices[1] = 0
        elif self.indices[0] < 1:
            self.indices[0] = 1
        for mrna in self.mrnas:
            mrna.clean_up_indices()

    def remove_invalid_features(self):
        empty_mrnas = []
        for i in xrange(len(self.mrnas)):
            if self.mrnas[i].indices[0] == 0:
                empty_mrnas.append(i)
        for j in reversed(empty_mrnas):
            self.mrnas.pop(j)
        for mrna in self.mrnas:
            mrna.remove_invalid_features()
        # TODO what if mrna.remove_empty leaves mrna with
        # no features? do we want to trash it?
            

    def length_of_shortest_cds_segment(self):
        min_length = self.mrnas[0].length_of_shortest_cds_segment()
        if len(self.mrnas) == 1:
            return min_length
        else:
            for mrna in self.mrnas:
                if mrna.length_of_shortest_cds_segment() < min_length:
                    min_length = mrna.length_of_shortest_cds_segment()
        return min_length

    def collidesRange(self, start, stop):
        if start <= self.indices[1] and stop >= self.indices[0]:
            return True
        return False

    def adjust_indices(self, n):
        self.indices = [i + n for i in self.indices]
        for mrna in self.mrnas:
            mrna.adjust_indices(n) 

    def trim(self, new_indices):
        self.trim_end(new_indices[1])
        self.trim_begin(new_indices[0])
        for mrna in self.mrnas:
            mrna.adjust_phase()
        self.clean_up_indices()
        self.remove_invalid_features()

    def to_gff(self):
        result = self.seq_name + "\t" + self.source + "\t"
        result += 'gene' + "\t" + str(self.indices[0]) + "\t"
        result += str(self.indices[1]) + "\t" + self.get_score()
        result += "\t" + self.strand + "\t" + "." + "\t"
        result += "ID=" + str(self.identifier) + ";Name=" + self.name + "\n"
        for mrna in self.mrnas:
            result += mrna.to_gff(self.seq_name, self.source, self.strand)
        return result




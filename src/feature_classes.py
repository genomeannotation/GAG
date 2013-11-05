#!/usr/bin/env python

import math
import sys

def length_of_segment(index_pair):
    return math.fabs(index_pair[1] - index_pair[0]) + 1

def adjust_index_pair(index_pair, n):
    return [i + n for i in index_pair]

class GenePart:
    def __init__(self, feature_type=None, identifier=None, name=None, indices=None, score=None, parent_id=None):
        self.feature_type = feature_type
        self.identifier = []
        if identifier is not None:
            self.identifier.append(identifier)
        self.name = []
        if name is not None:
            self.name.append(name)
        self.indices = []
        if indices is not None:
            self.indices.append(indices)
        self.score = []
        if score is not None:
            self.score.append(score)
        self.parent_id = parent_id

    def add_indices(self, ind):
        if isinstance(ind, list) and len(ind) is 2:
            self.indices.append(ind)
        else:
            raise ValueError()

    def add_name(self, name):
        self.name.append(name)

    def add_identifier(self, identifier):
        self.identifier.append(identifier)

    def add_score(self, score):
        self.score.append(score)

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

    def collidesRange(self, start, stop):
        for index_pair in self.indices:
            if start <= index_pair[1] and stop >= index_pair[0]:
                return True
        return False

    def remove_segment(self, segindex):
        try:
            self.identifier.pop(segindex)
            self.name.pop(segindex)
            self.indices.pop(segindex)
        except IndexError:
            sys.stderr.write("Trying to remove nonexistent segment " + str(segindex) + " from " + str(self))
        if len(self.score) > segindex:
            self.score.pop(segindex)

    def adjust_indices(self, n):
        self.indices = [adjust_index_pair(pair, n) for pair in self.indices]

    def trim_end(self, endindex):
        segs_to_cut = []
        for i, index_pair in enumerate(self.indices):
            if index_pair[0] > endindex:
                segs_to_cut.append(i)
            elif index_pair[1] > endindex:
                self.indices[i][1] = endindex
        for i in reversed(segs_to_cut):
            self.remove_segment(i)

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
        if len(self.identifier) <= i or self.parent_id is None:
            return None
        entry = "ID=" + str(self.identifier[i]) + ";"
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

    def __init__(self, identifier=None, name=None, indices=None, score=None, phase=None, parent_id=None):
        GenePart.__init__(self, feature_type='CDS', identifier=identifier, name=name, indices=indices, score=score, parent_id=parent_id)
        self.phase = []
        if phase is not None:
            self.phase.append(phase)
        self.annotations = []

    def get_phase(self, i):
        if self.phase and len(self.phase) > i:
            return self.phase[i]
        else:
            return "."

    def add_phase(self, ph):
        self.phase.append(ph)


class Exon(GenePart):

    def __init__(self, **kwargs):
        kwargs['feature_type'] = 'exon'
        GenePart.__init__(self, **kwargs)
        self.annotations = []


class MRNA:

    def __init__(self, identifier, name, indices, parent_id):
        self.identifier = identifier
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
        if self.exon:
            self.exon.adjust_indices(n)
        if self.cds:
            self.cds.adjust_indices(n)
        for feature in self.other_features:
            feature.adjust_indices(n)

    def trim_end(self, endindex):
        if self.indices[0] > endindex:
            print "foo"
            self.identifier = None
            self.name = None
            self.indices = None
            self.parent_id = None
            self.exon = None
            self.cds = None
            self.other_features = None
            # TODO write RIP message to stderr?
        elif self.indices[1] > endindex:
            self.indices[1] = endindex
            if self.cds:
                self.cds.trim_end(endindex)
            if self.exon:
                self.exon.trim_end(endindex)
            for feat in self.other_features:
                feat.trim_end(endindex)

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
            self.seq_name = None
            self.source = None
            self.indices = None 
            self.score = None
            self.strand = None 
            self.identifier = None 
            self.name = None 
            self.mrnas = None
        elif self.indices[1] > endindex:
            self.indices[1] = endindex
            for mrna in self.mrnas:
                mrna.trim_end(endindex)

    # beginindex is the new start index of sequence
    def trim_begin(self, beginindex):
        self.adjust_indices(-beginindex + 1)
            

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
        if n < 0 and math.fabs(n) > self.indices[0]:
            raise IndexError()
        else:
            self.indices = [i + n for i in self.indices]
            for mrna in self.mrnas:
                mrna.adjust_indices(n) 

    def to_gff(self):
        result = self.seq_name + "\t" + self.source + "\t"
        result += 'gene' + "\t" + str(self.indices[0]) + "\t"
        result += str(self.indices[1]) + "\t" + self.get_score()
        result += "\t" + self.strand + "\t" + "." + "\t"
        result += "ID=" + str(self.identifier) + ";Name=" + self.name + "\n"
        for mrna in self.mrnas:
            result += mrna.to_gff(self.seq_name, self.source, self.strand)
        return result




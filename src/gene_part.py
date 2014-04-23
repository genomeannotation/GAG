#!/usr/bin/env python

import math
import src.translator as translate

def length_of_segment(index_pair):
    return math.fabs(index_pair[1] - index_pair[0]) + 1

def adjust_index_pair(index_pair, n):
    return [i + n for i in index_pair]

def get_reversed_indices(indices):
    indices.reverse()
    [ind.reverse() for ind in indices]
    return indices

def one_line_indices_entry(indices, has_start, has_stop, is_cds):
    output = ""
    if not has_start:
        output += "<"
    output += str(indices[0]) + "\t"
    if not has_stop:
        output += ">"
    output += str(indices[1]) + "\t"
    if is_cds:
        output += "CDS\n"
    else:
        output += "mRNA\n"
    return output

def write_tbl_entry(indices, strand, has_start, has_stop, is_cds, annotations = None, phase=0):
    if not annotations:
        annotations = []
    output = ""
    if strand == "-":
        indices = get_reversed_indices(indices)
    if len(indices) == 1:
        output += one_line_indices_entry(indices[0], has_start, has_stop, is_cds)
    else:
        # Write first line of coordinates
        if not has_start:
            output += "<"
        output += str(indices[0][0]) + "\t" + str(indices[0][1]) + "\t" 
        if is_cds:
            output += "CDS\n"
        else:
            output += "mRNA\n"
        # Write middle lines
        for index_pair in indices[1:-1]:
            output += str(index_pair[0]) + "\t" + str(index_pair[1]) + "\n"
        # Write last line of coordinates
        output += str(indices[-1][0]) + "\t"
        if not has_stop:
            output += ">"
        output += str(indices[-1][1]) + "\n"
    if is_cds:
        output += "\t\t\tcodon_start\t" + str(phase+1) + "\n"
    output += "\t\t\tproduct\thypothetical protein\n"
    
    # Write the annotations
    for annot in annotations:
        output += '\t\t\t'+annot[0]+'\t'+annot[1]+'\n'
    
    return output

class GenePart:
    def __init__(self, feature_type=None, identifier=None,\
                 indices=None, score=None, strand='+', parent_id=None):
        self.feature_type = feature_type
        self.identifier = []
        if identifier is not None:
            self.identifier.append(identifier)
        self.indices = []
        if indices is not None:
            self.indices.append(indices)
        self.score = []
        if score is not None:
            self.score.append(score)
        self.strand = strand # Defauts to positive strand (?)
        self.parent_id = parent_id
        self.annotations = []

    def __str__(self):
        result = self.feature_type + " (first ID="
        result += str(self.identifier[0])+ ")"
        return result

    def add_indices(self, ind):
        if isinstance(ind, list) and len(ind) is 2:
            self.indices.append(ind)
            self.indices.sort()
            # TODO wait don't we need to sync this with the phase list for CDS?
        else:
            raise ValueError()

    def add_identifier(self, identifier):
        self.identifier.append(identifier)

    def add_score(self, score):
        self.score.append(score)
        
    def add_annotation(self, key, value):
        self.annotations.append([key, value])

    def gagflagged(self):
        for anno in self.annotations:
            if anno[0] == "gag_flag":
                return True
        return False

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

    def adjust_indices(self, n, start_index=1):
        for i, index_pair in enumerate(self.indices):
            if index_pair[0] >= start_index:
                self.indices[i] = adjust_index_pair(self.indices[i], n)
            elif index_pair[1] >= start_index:
                self.indices[i][1] += n

    def generate_attribute_entry(self, i):
        if len(self.identifier) <= i or self.parent_id is None:
            return None
        entry = "ID=" + str(self.identifier[i]) + ";"
        entry += "Parent=" + str(self.parent_id)
        for annot in self.annotations:
            entry += ';'+annot[0]+'='+annot[1]
        entry += '\n'
        return entry

    def to_gff(self, seq_name, source):
        result = ""
        for i in range(len(self.indices)):
            result += seq_name + "\t" + source + "\t"
            result += self.feature_type + "\t" + str(self.indices[i][0])
            result += "\t" + str(self.indices[i][1]) + "\t"
            result += str(self.get_score(i)) + "\t" + self.strand + "\t"
            result += str(self.get_phase(i)) + "\t"
            result += self.generate_attribute_entry(i)
        return result


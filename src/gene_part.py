#!/usr/bin/env python

import math
import copy
import sys
import src.translate as translate

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
                 indices=None, score=None, parent_id=None):
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
        self.parent_id = parent_id
        self.annotations = []

    def __str__(self):
        result = self.feature_type + " (first ID="
        result += str(self.identifier[0])+ ")"
        return result

    def add_indices(self, ind):
        if isinstance(ind, list) and len(ind) is 2:
            self.indices.append(ind)
        else:
            raise ValueError()

    def add_identifier(self, identifier):
        self.identifier.append(identifier)

    def add_score(self, score):
        self.score.append(score)
        
    def add_annotation(self, key, value):
        self.annotations.append([key, value])

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

    def remove_segment(self, segindex):
        try:
            self.identifier.pop(segindex)
            self.indices.pop(segindex)
        except IndexError:
            sys.stderr.write("Trying to remove nonexistent segment " + \
                             str(segindex) + " from " + str(self))
        if len(self.score) > segindex:
            self.score.pop(segindex)
        if self.feature_type == 'CDS' and len(self.phase) > segindex:
            self.phase.pop(segindex)

    def remove_trimmed_segments(self):
        segs_to_trim = []
        for i in xrange(len(self.indices)):
            if self.indices[i][0] == 0:
                segs_to_trim.append(i)
        for j in reversed(segs_to_trim):
            self.remove_segment(j)

    def adjust_indices(self, n, start_index=1):
        for i, index_pair in enumerate(self.indices):
            if index_pair[0] >= start_index:
                self.indices[i] = adjust_index_pair(self.indices[i], n)
            elif index_pair[1] >= start_index:
                self.indices[i][1] += n

    # Returns true if it survives invalidation. Returns false if the mRNA it lives on needs to die
    def invalidate_region(self, start, stop):
        for index in self.indices:
            # Index range contained in invalid region, 
            # mark for removal
            if start <= index[0] and stop >= index[1]:
                return False
            # Invalid region is in the middle of the index range, 
            # mark for removal
            elif start > index[0] and stop < index[1]:
                return False
            # The beginning is in the invalid region, 
            # trim beginning forward to invalid sequence stop
            elif start <= index[0] and stop >= index[0]:
                index[0] = stop+1
            # The end is in the invalid region, 
            # trim end back to invalid seq start
            elif start <= index[1] and stop >= index[1]:
                index[1] = start-1
        return True

    def valid_codon(self):  
        if length_of_segment(self.indices[0]) == 3:
            return True
        else:
            return False

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
        entry += "Parent=" + str(self.parent_id)
        for annot in self.annotations:
            entry += ';'+annot[0]+'='+annot[1]
        entry += '\n'
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

    def __init__(self, identifier=None, indices=None, \
                 score=None, phase=None, parent_id=None):
        GenePart.__init__(self, feature_type='CDS', identifier=identifier, \
                indices=indices, score=score, parent_id=parent_id)
        self.phase = []
        if phase is not None:
            self.phase.append(phase)

    def get_phase(self, i):
        if self.phase and len(self.phase) > i:
            return self.phase[i]
        else:
            return "."

    def add_phase(self, ph):
        self.phase.append(ph)

    def adjust_phase(self):
        for i in range(len(self.phase)):
            if self.indices[i][0] < 1:
                self.phase[i] = (self.phase[i] + self.indices[i][0] + -1) %3

    # Returns first and third indices regardless of whether 
    # CDS actually has a start codon
    def get_start_indices(self, phase):
        if phase == '+':
            first_index = self.indices[0][0]
            return [first_index, first_index+2]
        elif phase == '-':
            first_index = self.indices[0][1]
            return [first_index-2, first_index]

    def get_stop_indices(self, phase):
        if phase == '+':
            last_index_pair = self.indices[len(self.indices)-1]
            last_index = last_index_pair[1]
            return [last_index-2, last_index]
        elif phase == '-':
            last_index_pair = self.indices[len(self.indices)-1]
            last_index = last_index_pair[0]
            return [last_index, last_index+2]

    # Override to account for phase
    def invalidate_region(self, start, stop):
        for i, index in enumerate(self.indices):
            # index range contained in invalid region, 
            # mark for removal
            if start <= index[0] and stop >= index[1]:
                return False
            # Invalid region is in the middle of the index range, 
            # mark for removal
            elif start > index[0] and stop < index[1]:
                return False
            # The beginning is in the invalid region, 
            # trim beginning forward to invalid sequence stop
            elif start <= index[0] and stop >= index[0]:
                self.phase[0] = (self.phase[0] - ((stop+1)-index[0])%3)%3 # adjust phase
                index[0] = stop+1
            # The end is in the invalid region, 
            # trim end back to invalid seq start
            elif start <= index[1] and stop >= index[1]:
                index[1] = start-1
        return True

    def extract_sequence(self, seq_object, strand):
        seq = ''
        if strand == '+':
            for i in xrange(len(self.indices)):
                index_pair = self.indices[i]
                subseq = seq_object.get_subseq(index_pair[0], index_pair[1])
                if subseq:
                    seq += subseq
        elif strand == '-':
            for i in xrange(len(self.indices)):
                index_pair = self.indices[i]
                non_reversed_seq = seq_object.get_subseq(index_pair[0], index_pair[1])
                if non_reversed_seq:
                    seq += translate.reverse_complement(non_reversed_seq)
        return seq

    def to_tbl(self, strand, has_start, has_stop):
        indices = copy.deepcopy(self.indices)
        phase = self.phase[0]
        return write_tbl_entry(indices, strand, has_start, has_stop, True, self.annotations, phase)


class Exon(GenePart):

    def __init__(self, **kwargs):
        kwargs['feature_type'] = 'exon'
        GenePart.__init__(self, **kwargs)

    def to_tbl(self, strand, has_start, has_stop):
        indices = copy.deepcopy(self.indices)
        return write_tbl_entry(indices, strand, has_start, has_stop, False, self.annotations)


#!/usr/bin/env python

import copy
from src.gene_part import *

class CDS(GenePart):

    def __init__(self, identifier=None, indices=None, \
                 score=None, phase=None, strand=None, parent_id=None):
        GenePart.__init__(self, feature_type='CDS', identifier=identifier, \
                indices=indices, score=score, strand=strand, parent_id=parent_id)
        self.phase = []
        if phase is not None:
            self.phase.append(phase)

    def get_phase(self, i):
        """Returns phase for given segment of CDS."""
        if self.phase and len(self.phase) > i:
            return self.phase[i]
        else:
            return "."

    def add_phase(self, ph):
        """Appends phase to CDS"""
        self.phase.append(ph)

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

    def to_tbl(self, has_start, has_stop):
        indices = copy.deepcopy(self.indices)
        phase = self.phase[0]
        return write_tbl_entry(indices, self.strand, has_start, has_stop, True, self.annotations, phase)


#!/usr/bin/env python

import math
import sys
from feature_tbl_entry import FeatureTblEntry

def length_of_segment(index_pair):
    return math.fabs(index_pair[1] - index_pair[0]) + 1

def trimmed_completely(gene_inds, seq_inds):
    return seq_inds == [0,0] or gene_inds[0] > seq_inds[1] or gene_inds[1] < seq_inds[0]

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

    def is_empty(self):
        return self.indices == [0, 0]

    def length(self):
        return length_of_segment(self.indices)

    def get_score(self):
        if self.score:
            return self.score
        else:
            return '.'

    def add_mrna(self, mrna):
        self.mrnas.append(mrna)

    def contains_mrna_named(self, name):
        for mrna in self.mrnas:
            if mrna.name == name:
                return True
        return False

    def remove_mrnas_with_cds_shorter_than(self, min_length):
        # TODO for now this also removes mrnas with NO cds, but do we want that?
        to_remove = []
        if self.mrnas:
            for mrna in self.mrnas:
                if mrna.cds:
                    if mrna.cds.length() < min_length:
                        to_remove.append(mrna)
                else:
                    to_remove.append(mrna)
        for m in to_remove:
            self.mrnas.remove(m)

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
        self.adjust_indices(-beginindex + 1, 1)

    def clean_up_indices(self):
        if self.indices[1] < 1:
            self.indices[0] = 0
            self.indices[1] = 0
        elif self.indices[0] < 1:
            self.indices[0] = 1
        for mrna in self.mrnas:
            mrna.clean_up_indices()

    def create_starts_and_stops(self, fasta):
        for mrna in self.mrnas:
            mrna.create_start_and_stop_if_necessary(fasta, self.seq_name, self.strand)

    def remove_first_cds_segment_if_shorter_than(self, min_length):
        if self.mrnas:
            for mrna in self.mrnas:
                mrna.remove_first_cds_segment_if_shorter_than(min_length)

    def remove_invalid_features(self):
        # remove mrnas with indices[0] == 0
        self.mrnas = [m for m in self.mrnas if m.indices[0] != 0]
        for mrna in self.mrnas:
            mrna.remove_invalid_features()

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

    def adjust_indices(self, n, start_index=1):
        if self.indices[0] >= start_index:
            self.indices = [i + n for i in self.indices]
        elif self.indices[1] >= start_index:
            self.indices[1] += n
        for mrna in self.mrnas:
            mrna.adjust_indices(n, start_index) 

    def trim(self, new_indices):
        if trimmed_completely(self.indices, new_indices):
            self.mrnas = []
            self.indices = [0, 0]
        else:
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

    def to_tbl_entries(self, annotator):
        entries = []
        geneEntry = FeatureTblEntry()
        geneEntry.set_type("gene")
        geneEntry.set_name(self.name)
        geneEntry.set_seq_name(self.seq_name)
        geneEntry.add_coordinates(self.indices[0], self.indices[1])
        geneEntry.set_strand(self.strand)
        geneEntry.set_phase(0)
        geneEntry.set_partial_info(True, True) # Pretend there's a start and stop codon for genes
        annotator.annotate_gene(geneEntry)
        entries.append(geneEntry)

        hypothetical = False
        for mrna in self.mrnas: 
            mrna_entries = mrna.to_tbl_entries(annotator, self.strand)
            for mrna_entry in mrna_entries:
                mrna_entry.set_seq_name(self.seq_name)
                if mrna_entry.is_hypothetical():
                    hypothetical = True
                entries.append(mrna_entry)

        # Hypothetical genes don't get gene names
        if hypothetical == True:
            geneAnnot = geneEntry.get_annotation('gene')
            geneEntry.remove_annotation('gene')
            if geneAnnot:
                geneEntry.add_annotation('note', 'gene '+geneAnnot)

        return entries




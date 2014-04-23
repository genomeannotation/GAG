#!/usr/bin/env python

import math

def length_of_segment(index_pair):
    return math.fabs(index_pair[1] - index_pair[0]) + 1

class Gene:

    def __init__(self, seq_name, source, indices, strand, identifier, annotations=None, score=None):
        self.seq_name = seq_name
        self.source = source
        self.indices = indices
        self.score = score
        self.strand = strand
        self.identifier = identifier
        self.mrnas = []
        if not annotations:
            self.annotations = []
        else:
            self.annotations = annotations
        self.death_flagged = False

    def __str__(self):
        result = "Gene (ID=" + str(self.identifier) 
        result += ", seq_name=" + self.seq_name
        result += ") containing " + str(len(self.mrnas))
        result += " mrnas"
        return result

#    def get_valid_mrnas(self):
#        return [mrna for mrna in self.mrnas if not mrna.death_flagged]

    def get_mrna_ids(self):
        result = []
        for mrna in self.mrnas:
            result.append(mrna.identifier)
        return result
        
    def add_annotation(self, key, value):
        self.annotations.append([key, value])
        
    def length(self):
        return length_of_segment(self.indices)

    def gagflagged(self):
        for anno in self.annotations:
            if anno[0] == "gag_flag":
                return True
        return False

    def number_of_gagflags(self):
        total = 0
        for mrna in self.mrnas:
            total += mrna.number_of_gagflags()
        if self.gagflagged():
            total += 1
        return total

    def get_score(self):
        if self.score:
            return self.score
        else:
            return '.'

    def get_longest_exon(self):
        longest = 0
        for mrna in self.mrnas:
            length = mrna.get_longest_exon()
            if length > longest:
                longest = length
        return longest

    def get_shortest_exon(self):
        shortest = 0
        for mrna in self.mrnas:
            length = mrna.get_shortest_exon()
            if shortest == 0 or length < shortest:
                shortest = length
        return shortest

    def get_total_exon_length(self):
        total = 0
        for mrna in self.mrnas:
            total += mrna.get_total_exon_length()
        return total
    
    def get_num_exons(self):
        total = 0
        for mrna in self.mrnas:
            total += mrna.get_num_exons()
        return total

    def get_longest_intron(self):
        longest = 0
        for mrna in self.mrnas:
            length = mrna.get_longest_intron()
            if length > longest:
                longest = length
        return longest

    def get_shortest_intron(self):
        shortest = 0
        for mrna in self.mrnas:
            length = mrna.get_shortest_intron()
            if shortest == 0 or length < shortest:
                shortest = length
        return shortest

    def get_total_intron_length(self):
        total = 0
        for mrna in self.mrnas:
            total += mrna.get_total_intron_length()
        return total

    def get_num_introns(self):
        total = 0
        for mrna in self.mrnas:
            total += mrna.get_num_introns()
        return total
    
    def create_starts_and_stops(self, seq_object):
        for mrna in self.mrnas:
            mrna.create_start_and_stop_if_necessary(seq_object, self.strand)

    def adjust_indices(self, n, start_index=1):
        """Adds 'n' to both indices, checking to ensure that they fall after an optional start index"""
        if self.indices[0] >= start_index:
            self.indices = [i + n for i in self.indices]
        elif self.indices[1] >= start_index:
            self.indices[1] += n
        for mrna in self.mrnas:
            mrna.adjust_indices(n, start_index)

    def get_partial_info(self):
        results = {"complete": 0, "start_no_stop": 0, "stop_no_start": 0, "no_stop_no_start": 0}
        for mrna in self.mrnas:
            if mrna.has_start():
                if mrna.has_stop():
                    results["complete"] += 1
                else:
                    results["start_no_stop"] += 1
            else:
                # No start ...
                if mrna.has_stop():
                    results["stop_no_start"] += 1
                else:
                    results["no_stop_no_start"] += 1
        return results

    def remove_mrnas_with_internal_stops(self, seq_helper):
        for mrna in self.mrnas:
            if seq_helper.mrna_contains_internal_stop(mrna):
                mrna.death_flagged = True
        self.mrnas = [m for m in self.mrnas if not m.death_flagged]

    def contains_mrna(self, mrna_id):
        for mrna in self.mrnas:
            if mrna.identifier == mrna_id:
                return True
        return False

    def cds_to_gff(self, seq_id, mrna_id):
        for mrna in self.mrnas:
            if mrna.identifier == mrna_id and mrna.cds:
                return mrna.cds_to_gff(seq_id, self.source)
        return ""

    def cds_to_tbl(self, mrna_id):
        for mrna in self.mrnas:
            if mrna.identifier == mrna_id and mrna.cds:
                return mrna.cds_to_tbl()
        return ""

    def to_mrna_fasta(self, seq_helper):
        result = ""
        for mrna in self.mrnas:
            result += seq_helper.mrna_to_fasta(mrna)
        return result

    def to_cds_fasta(self, seq_helper):
        result = ""
        for mrna in self.mrnas:
            result += seq_helper.mrna_to_cds_fasta(mrna)
        return result

    def to_protein_fasta(self, seq_helper):
        result = ""
        for mrna in self.mrnas:
            result += seq_helper.mrna_to_protein_fasta(mrna)
        return result

    def to_gff(self):
        result = self.seq_name + "\t" + self.source + "\t"
        result += 'gene' + "\t" + str(self.indices[0]) + "\t"
        result += str(self.indices[1]) + "\t" + self.get_score()
        result += "\t" + self.strand + "\t" + "." + "\t"
        result += "ID=" + str(self.identifier)
        for annot in self.annotations:
            result += ';'+annot[0]+'='+annot[1]
        result += '\n'
        for mrna in self.mrnas:
            result += mrna.to_gff(self.seq_name, self.source)
        return result

    def to_tbl(self):
        if self.strand == "-":
            indices = [self.indices[1], self.indices[0]]
        else:
            indices = self.indices
        output = str(indices[0]) + "\t" + str(indices[1]) + "\t" + "gene\n"
        output += "\t\t\tlocus_tag\t" + self.identifier + "\n"
        for annot in self.annotations:
            output += '\t\t\t'+annot[0]+'\t'+annot[1]+'\n'
        for mrna in self.mrnas:
            output += mrna.to_tbl()
        return output




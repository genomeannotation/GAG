#!/usr/bin/env python


class SeqHelper:

    def __init__(self, bases):
        self.full_sequence = bases

    def mrna_to_fasta(self, mrna):
        """Writes a two-line fasta-style entry consisting of all exonic sequence."""

        identifier = mrna.identifier
        strand = mrna.strand
        indices = mrna.exon.indices
        return self.list_of_index_pairs_to_fasta(identifier, strand, indices)

    def mrna_to_cds_fasta(self, mrna):
        """Writes a two-line fasta-style entry consisting of all CDS sequence."""

        identifier = mrna.identifier + " CDS"
        strand = mrna.strand
        indices = mrna.cds.indices
        return self.list_of_index_pairs_to_fasta(identifier, strand, indices)

    def list_of_index_pairs_to_fasta(self, identifier, strand, indices):
        if strand == '+':
            positive = True
        else:
            positive = False
        result = identifier + "\n"
        for index_pair in indices:
            start = index_pair[0]-1
            stop = index_pair[1]
            if positive:
                result += self.full_sequence[start:stop]
            else:
                result += self.full_sequence[start:stop][::-1]
        return result + "\n"


#!/usr/bin/env python

from src.translate import translate

class SeqHelper:

    def __init__(self, bases):
        self.full_sequence = bases

    def mrna_to_fasta(self, mrna):
        """Writes a two-line fasta-style entry consisting of all exonic sequence."""

        identifier = mrna.identifier
        strand = mrna.strand
        indices = mrna.exon.indices
        return self.id_and_indices_to_fasta(identifier, strand, indices)

    def mrna_to_cds_fasta(self, mrna):
        """Writes a two-line fasta-style entry consisting of all CDS sequence."""

        identifier = mrna.identifier + " CDS"
        strand = mrna.strand
        indices = mrna.cds.indices
        return self.id_and_indices_to_fasta(identifier, strand, indices)

    def mrna_to_protein_fasta(self, mrna):
        """Writes a two-line fasta-style entry consisting of the translation of CDS sequence."""

        identifier = mrna.identifier + " protein"
        strand = mrna.strand
        indices = mrna.cds.indices
        untranslated = self.get_sequence_from_indices(strand, indices)
        return identifier + "\n" + translate(untranslated, strand) + "\n"

    def id_and_indices_to_fasta(self, identifier, strand, indices):
        result = identifier + "\n"
        result += self.get_sequence_from_indices(strand, indices) + "\n"
        return result

    def get_sequence_from_indices(self, strand, indices):
        if strand == '+':
            positive = True
        else:
            positive = False
        result = ""
        for index_pair in indices:
            start = index_pair[0]-1
            stop = index_pair[1]
            if positive:
                result += self.full_sequence[start:stop]
            else:
                result += self.full_sequence[start:stop][::-1]
        return result

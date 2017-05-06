#!/usr/bin/env python
# coding=utf-8

from src.translator import translate, reverse_complement, contains_internal_stop


class SeqHelper(object):
    def __init__(self, bases):
        self.full_sequence = bases

    def mrna_contains_internal_stop(self, mrna):
        if not mrna.cds:
            return False
        strand = mrna.strand
        indices = mrna.cds.indices
        sequence = self.get_sequence_from_indices(strand, indices)
        return contains_internal_stop(sequence, strand)

    def mrna_to_fasta(self, mrna):
        """Writes a two-line fasta-style entry consisting of all exonic sequence."""
        if not mrna.exon:
            return ""
        identifier = ">" + mrna.identifier
        metadata = "ID=" + mrna.identifier + "|Parent=" + mrna.parent_id + "|Name=" + mrna.name
        strand = mrna.strand
        indices = mrna.exon.indices
        return self.id_and_indices_to_fasta(identifier, strand, indices, metadata)

    def mrna_to_cds_fasta(self, mrna):
        """Writes a two-line fasta-style entry consisting of all CDS sequence."""
        if not mrna.cds:
            return ""
        identifier = ">CDS|" + mrna.identifier
        metadata = "ID=" + mrna.identifier + "|Parent=" + mrna.parent_id + "|Name=" + mrna.name
        strand = mrna.strand
        indices = mrna.cds.indices
        return self.id_and_indices_to_fasta(identifier, strand, indices, metadata)

    def mrna_to_protein_fasta(self, mrna):
        """Writes a two-line fasta-style entry consisting of the translation of CDS sequence."""
        if not mrna.cds:
            return ""
        identifier = ">protein|" + mrna.identifier
        metadata = "ID=" + mrna.identifier + "|Parent=" + mrna.parent_id + "|Name=" + mrna.name
        strand = mrna.strand
        indices = mrna.cds.indices
        untranslated = self.get_sequence_from_indices(strand, indices)
        # Take phase into account
        if strand is "+":
            ph = mrna.cds.get_phase(0)
        elif strand is "-":
            ph = mrna.cds.get_phase(-1)
        else:
            ph = 0
        untranslated = untranslated[ph:]
        return identifier + " " + metadata + "\n" + translate(untranslated, "+") + "\n"

    def id_and_indices_to_fasta(self, identifier, strand, indices, metadata=None):
        result = identifier
        if metadata:
            result += " " + metadata + "\n"
        else:
            result += "\n"
        result += self.get_sequence_from_indices(strand, indices) + "\n"
        return result

    def get_sequence_from_indices(self, strand, indices):
        result = ""
        for index_pair in indices:
            start = index_pair[0] - 1
            stop = index_pair[1]
            result += self.full_sequence[start:stop]
        if strand == '-':
            result = reverse_complement(result)
        return result

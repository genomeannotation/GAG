#!/usr/bin/env python


class SeqHelper:

    def __init__(self, bases):
        self.full_sequence = bases

    def mrna_to_fasta(self, mrna):
        """Writes a two-line fasta-style entry consisting of all exonic sequence."""

        if mrna.strand == '+':
            positive = True
        else:
            positive = False
        result = mrna.identifier + "\n"
        for index_pair in mrna.exon.indices:
            start = index_pair[0]-1
            stop = index_pair[1]
            if positive:
                result += self.full_sequence[start:stop]
            else:
                result += self.full_sequence[start:stop][::-1]
        return result + "\n"

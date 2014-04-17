#!/usr/bin/env python


class SeqHelper:

    def __init__(self, bases):
        self.full_sequence = bases

    def mrna_to_fasta(self, mrna):
        """Writes a two-line fasta-style entry consisting of all exonic sequence."""

        result = mrna.identifier + "\n"
        for index_pair in mrna.exon.indices:
            start = index_pair[0] - 1
            stop = index_pair[1]
            result += self.full_sequence[start:stop]
        return result + "\n"

#!/usr/bin/env python

class Sequence:

    def __init__(self, header="", bases=""):
        self.header = header
        self.bases = bases
        self.genes = []

    def __str__(self):
        result = "Sequence " + self.header
        result += " of length " + str(len(self.bases))
        result += " containing "
        result += str(len(self.genes))
        result += " genes\n"
        return result

    def add_gene(self, gene):
        self.genes.append(gene)

    def remove_gene(self, gene_id):
        self.genes = [g for g in self.genes if g.identifier != gene_id]

    def to_fasta(self):
        result = '>' + self.header + '\n'
        result += self.bases + '\n'
        return result

    # Given a position in the sequence, returns the number of Ns 
    # from that position forward 
    # (returns 0 if the base at that position is not N)
    def how_many_Ns_forward(self, position):
        index = position-1
        if self.bases[index] != 'N' and self.bases[index] != 'n':
            return 0
        else:
            count = 1
            index += 1
            for base in self.bases[index:]:
                if base != 'N' and base != 'n':
                    return count
                else:
                    count += 1
            return count

    # Given a position in the fasta, returns the number of Ns 
    # from that position backward 
    # (returns 0 if the base at that position is not N)
    def how_many_Ns_backward(self, position):
        index = position-1
        if self.bases[index] != 'N' and self.bases[index] != 'n':
            return 0
        else:
            count = 1
            index -= 1
            for base in self.bases[index::-1]:
                if base != 'N' and base != 'n':
                    return count
                else:
                    count += 1
            return count

    def trim_region(self, start, stop):
        if stop > len(self.bases):
            sys.stderr.write("Sequence.trim called on sequence that is too short;"+\
                    " doing nothing.")
            return
        self.bases = self.bases[:start-1] + self.bases[stop:]
        offset = -(stop - start + 1)
        to_remove = []
        # TODO don't autoremove them; adjust indices and verify
        for gene in self.genes:
            for index in gene.indices:
                if index >= start and index <= stop:
                    to_remove.append(gene)
        self.genes = [g for g in self.genes if g not in to_remove]
        
    def get_subseq(self, start, stop):
        if stop > len(self.bases):
            return ""
        return self.bases[start-1:stop]

    def to_tbl(self):
        result = ""
        for gene in self.genes:
            result += gene.to_tbl()
        return result

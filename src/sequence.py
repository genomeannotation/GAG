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

    def contains_gene(self, gene_id):
        for gene in self.genes:
            if gene.identifier == gene_id:
                return True
        return False

    def contains_mrna(self, mrna_id):
        for gene in self.genes:
            for mrna in gene.mrnas:
                if mrna.identifier == mrna_id:
                    return True
        return False

    def get_locus_tag(self):
        for gene in self.genes:
            gene_id = str(gene.identifier)
            locus_tag = gene_id.split('_')[0]
            return locus_tag
        return ""

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

    def create_starts_and_stops(self):
        for gene in self.genes:
            gene.create_starts_and_stops(self)

    def remove_mrna(self, args):
        for gene in self.genes:
            gene.mrnas = [m for m in gene.mrnas if m.identifier not in args]

    def remove_gene(self, args):
        self.genes = [g for g in self.genes if g.identifier not in args]

    def remove_mrnas_with_cds_shorter_than(self, min_length):
        for gene in self.genes:
            gene.remove_mrnas_with_cds_shorter_than(min_length)

    def invalidate_region(self, start, stop):
        for gene in seq.genes:
            gene.invalidate_region(start, stop)

    def extract_cds_seq(self, mrna_id):
        for gene in self.genes:
            for mrna in gene.mrnas:
                if mrna.identifier == mrna_id and mrna.cds:
                    return mrna.cds.extract_sequence(self, gene.strand)

    def to_tbl(self):
        result = "1\t" + str(len(self.bases)) + "\tREFERENCE\n"
        result += "\t\t\tPBARC\t12345\n"
        for gene in self.genes:
            result += gene.to_tbl()
        return result

    def to_gff(self):
        result = ""
        for gene in self.genes:
            result += gene.to_gff()
        return result

    def gene_to_gff(self, gene_id):
        for gene in self.genes:
            if gene.identifier == gene_id:
                return gene.to_gff()
        return ""

    def gene_to_tbl(self, gene_id):
        for gene in self.genes:
            if gene.identifier == gene_id:
                return gene.to_tbl()
        return ""
        
    def stats(self):
        stats = dict()
        
        stats["seq_length"] = len(self.bases)
        
        return stats

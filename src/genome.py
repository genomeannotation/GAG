#!/usr/bin/env python

from feature_tbl_entry import FeatureTblEntry
from translate import *
import sys
import os.path

class Genome:

    def __init__(self):
        self.fasta = None
        self.gff = None
        self.annot = None
        self.entries = []

    def verify_file(self, filename):
        return os.path.exists(filename) 

    def addEntry(self, entry):
        entries.append(entry)
    
    def addEntries(self, entries):
        [self.entries.append(entry) for entry in entries]

    # maybe should be called 'create_start_codon_GenePart if sequence contains start codon'
    def verify_start_codon(self, mrna, seq_id):
        indices = mrna.get_cds_indices()
        seq = self.fasta.get_subseq(seq_id, indices[0])
        if has_start_codon(seq):
            mrna.add_start_codon(indices[0][0])

    def verify_stop_codon(self, mrna, seq_id):
        indices = mrna.get_cds_indices()
        last_pair = indices[len(indices)-1]
        seq = self.fasta.get_subseq(seq_id, last_pair)
        if has_stop_codon(seq):
            mrna.add_stop_codon(last_pair[1])

    def verify_all_starts_and_stops(self):
        for gene in self.gff.genes:
            for mrna in gene.mrnas:
                if not mrna.has_start_codon():
                    self.verify_start_codon(mrna, gene.seq_name)
                if not mrna.has_stop_codon():
                    self.verify_stop_codon(mrna, gene.seq_name)

    def generateEntries(self):
        for gene in self.gff.genes:
            newEntries = gene.to_tbl_entries()
            for entry in newEntries:
                if entry.type == 'gene':
                    self.annot.annotate_gene(entry)
                elif entry.type == 'CDS':
                    self.annot.annotate_cds(entry)
                elif entry.type == 'mRNA':
                    self.annot.annotate_mrna(entry)
                self.entries.append(entry)

    def write_string(self, genes = None, errors = None):
        output = ''

        if self.fasta == None or self.gff == None or self.annot == None:
            return output
       
        output += '>Feature SeqId\n'

        for seq in self.fasta.entries:
            if genes != None and not genes:
                return output

            entries = []
            for gene in self.gff.genes:
                if gene.seq_name != seq[0]:
                    continue

                if genes != None and gene.name not in genes:
                    continue

                if genes != None:
                    genes.remove(gene.name)

                newEntries = gene.to_tbl_entries()
                for entry in newEntries:
                    if entry.type == 'gene':
                        self.annot.annotate_gene(entry)
                    elif entry.type == 'CDS':
                        self.annot.annotate_cds(entry)
                    elif entry.type == 'mRNA':
                        self.annot.annotate_mrna(entry)
                    entries.append(entry)

            # If there are any entries, write this section of the tbl file
            if len(entries) > 0:
                output += '>Feature '+seq[0]+'\n'
                output += '1\t'+str(len(seq[1]))+'\tREFERENCE\n\t\t\tPBARC\t12345\n'

                for entry in entries:    
                    output += entry.write_to_string()+'\n'
        return output

    def remove_all_gene_segments(self, prefix):
        self.gff.remove_all_gene_segments(prefix)

    def obliterate_genes_related_to_mrnas(self, mrna_names):
        self.gff.obliterate_genes_related_to_mrnas(mrna_names)

    def rename_maker_mrnas(self):
        count = 1000000
        for gene in self.gff.genes:
            for mrna in gene.mrnas:
                if mrna.is_maker_mrna():
                    old_name = mrna.name
                    new_name = 'BDOR_' + str(count)
                    mrna.name = new_name
                    self.annot.rename_mrna(old_name, new_name)
                    count += 1


#!/usr/bin/env python

from feature_tbl_entry import FeatureTblEntry
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


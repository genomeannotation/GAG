#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import csv
import cmd
import readline
from src.fasta import Fasta
from src.gene import Gene
from src.mrna import MRNA
from src.gene_part import GenePart, CDS, Exon
from src.gff import GFF
from src.bed import Bed
from src.feature_tbl_entry import FeatureTblEntry
from src.annotator import Annotator

class GagCmd(cmd.Cmd):

    def __init__(self):
        cmd.Cmd.__init__(self)
        self.prompt = "GAG> "

        readline.set_history_length(1000)
        readline.read_history_file('.gaghistory')
                

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_readgff(self):
        print("Usage: readgff <file_name>\n")
        print("Read the gff file. Any unsaved changes")
        print("to the currently loaded gff will be lost.\n")

    def do_readgff(self, line):
        self.gff = GFF()
        with open(line, 'rb') as gfffile:
            gffreader = csv.reader(gfffile, delimiter='\t')
            self.gff.read_file(gffreader)

    def help_readfasta(self):
        print("Usage: readfasta <file_name>\n")
        print("Read the fasta file. Any unsaved changes")
        print("to the currently loaded fasta will be lost.\n")

    def do_readfasta(self, line):
        self.fasta = Fasta()
        self.fasta.read_file(line)

    def help_readtrinotate(self):
        print("Usage: readtrinotate <file_name>\n")
        print("Read annotations from a trinotate file.\n")

    def do_readtrinotate(self, line):
        self.annot = Annotator()
        self.annot.read_from_file(line)

    def help_applybed(self):
        print("Usage: applybed <file_name>\n")
        print("Applies a bed file to the data. This will")
        print("trim out any sequences that aren't in the")
        print("ranges of the bed file and automagically")
        print("update the GFF file accordingly.\n")

    def do_applybed(self, line):
        bed = Bed()
        with open(line, 'rb') as bedfile:
            bedreader = csv.reader(bedfile, delimiter='\t')
            bed.read_file(bedreader)
            self.fasta.apply_bed(bed)
            self.gff.apply_bed(bed)
            self.gff.remove_empty_genes()

    def help_writetbl(self):
        print("Usage: writetbl <file_name>\n")
        print("Write a sweet feature table to the specified file.\n")

    def do_writetbl(self, line):
        outFile = open(line, 'w')
        outFile.write('>Feature SeqId\n')

        for seq in self.fasta.entries:
            outFile.write('>Feature '+seq[0]+'\n')
            outFile.write('1\t'+str(len(seq[1]))+'\tREFERENCE\n\t\t\tPBARC\t12345\n')
            for gene in self.gff.genes:
                if gene.seq_name != seq[0]:
                    continue

                entries = gene.to_tbl_entries()
                for entry in entries:    
                    if entry.type == 'gene':
                        self.annot.annotate_gene(entry)
                    elif entry.type == 'CDS':
                        self.annot.annotate_cds(entry)
                    elif entry.type == 'mRNA':
                        self.annot.annotate_mrna(entry)
                    outFile.write(entry.write_to_string()+'\n')
        outFile.close()

    def do_barfsession(self, line):
        if len(line) == 0:
            print("Usage: barfsession <directory>\n")
            return

        os.system('mkdir '+line)
        
        # Write the gff
        with open(line+'/gag.gff', 'w') as gff:
            for gene in self.gff.genes:
                gff.write(gene.to_gff())
            gff.close()

        # Write the fasta
        with open(line+'/gag.fasta', 'w') as fasta:
            fasta.write(self.fasta.write_string())
            fasta.close()

        # Write the annotations
        self.annot.write_to_file(line+'/gag.trinotate')

    def do_loadsession(self, line):
        # Read the gff
        self.do_readgff(line+'/gag.gff')

        # Read the fasta
        self.do_readfasta(line+'/gag.fasta')

        # Read the annotations
        self.do_readtrinotate(line+'/gag.trinotate')

    def do_barfgenegff(self, line):
        for gene in self.gff.genes:
            if gene.name == line:
                print(gene.to_gff())

    def do_barfgenetbl(self, line):
        for gene in self.gff.genes:
            if gene.name == line:
                entries = gene.to_tbl_entries()
                for entry in entries:    
                    if entry.type == 'gene':
                        self.annot.annotate_gene(entry)
                    elif entry.type == 'CDS':
                        self.annot.annotate_cds(entry)
                    elif entry.type == 'mRNA':
                        self.annot.annotate_mrna(entry)
                    print(entry.write_to_string()+'\n')
                return

    def help_exit(self):
        print("Exit this console.\n")

    def do_exit(self, line):
        return True

if __name__ == '__main__':
    GagCmd().cmdloop("Welcome to the GAG console!")

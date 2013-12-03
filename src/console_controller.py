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
from src.genome import Genome
from src.annotator import Annotator
from src.translate import translate

class ConsoleController:

    def __init__(self):
        self.genome = Genome()

    def read_gff(self, line):
        self.genome.gff = GFF()
        with open(line, 'rb') as gfffile:
            gffreader = csv.reader(gfffile, delimiter='\t')
            self.genome.gff.read_file(gffreader)
        return self.genome.gff

    def read_fasta(self, line):
        self.genome.fasta = Fasta()
        self.genome.fasta.read_file(line)
        return self.genome.fasta

    def barf_gff(self, line):
        for gene in self.genome.gff.genes:
            if gene.name == line:
                print(gene.to_gff())

"""

    def help_readtrinotate(self):
        print("Usage: readtrinotate <file_name>\n")
        print("Read annotations from a trinotate file.\n")

    def do_readtrinotate(self, line):
        self.genome.annot = Annotator()
        self.genome.annot.read_from_file(line)

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
            self.genome.fasta.apply_bed(bed)
            self.genome.gff.apply_bed(bed)
            self.genome.gff.remove_empty_genes()

    def help_writetbl(self):
        print("Usage: writetbl <file_name>\n")
        print("Write a sweet feature table to the specified file.\n")

    def do_writetbl(self, line):
        with open(line, 'w') as outFile:
            self.genome.write_file(outFile)
            outFile.close()

    def do_barferrsubset(self, line):
        args = line.split()

        if len(args) < 2:
            print("Usage: barferrsubset <directory> <errorcode>\n")
            return

        outdir = args[0]
        err = args[1]

    def do_barfsession(self, line):
        if len(line) == 0:
            print("Usage: barfsession <directory>\n")
            return

        os.system('mkdir '+line)
        
        # Write the gff
        with open(line+'/gag.gff', 'w') as gff:
            for gene in self.genome.gff.genes:
                gff.write(gene.to_gff())
            gff.close()

        # Write the fasta
        with open(line+'/gag.fasta', 'w') as fasta:
            fasta.write(self.genome.fasta.write_string())
            fasta.close()

        # Write the annotations
        self.genome.annot.write_to_file(line+'/gag.trinotate')

    def do_loadsession(self, line):
        # Read the gff
        self.do_readgff(line+'/gag.gff')

        # Read the fasta
        self.do_readfasta(line+'/gag.fasta')

        # Read the annotations
        self.do_readtrinotate(line+'/gag.trinotate')


    def do_barfgenetbl(self, line):
        self.genome.write_file(sys.stdout, set(line.split()))

    def do_barfseq(self, line):
        args = line.split(' ')
        print(str(self.genome.fasta.get_subseq(args[0], [int(args[1]), int(args[2])]))+'\n')

    def do_ducttapeseqframes(self, line):
        for gene in self.genome.gff.genes:
            for mrna in gene.mrnas:
                if mrna.name == line:
                    seq = self.genome.fasta.get_subseq(gene.seq_name, mrna.cds.indices[0])
                    pseq1 = translate(seq, 1, '+')
                    pseq2 = translate(seq, 2, '+')
                    pseq3 = translate(seq, 3, '+')
                    nseq1 = translate(seq, 1, '-')
                    nseq2 = translate(seq, 2, '-')
                    nseq3 = translate(seq, 3, '-')

                    pepSeq = self.genome.annot.get_entry(line)[9]
                    if pepSeq.find(pseq1) != -1:
                        print('+1')
                    elif pepSeq.find(pseq2) != -1:
                        print('+2')
                    elif pepSeq.find(pseq3) != -1:
                        print('+3')
                    elif pepSeq.find(nseq1) != -1:
                        print('-1')
                    elif pepSeq.find(nseq2) != -1:
                        print('-2')
                    elif pepSeq.find(nseq3) != -1:
                        print('-3')

                    return

    def help_exit(self):
        print("Exit this console.\n")

    def do_exit(self, line):
        return True

if __name__ == '__main__':
    GagCmd().cmdloop("Welcome to the GAG console!")
"""

#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import csv
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

## Setup, loading and saving sessions

    def __init__(self):
        self.genome = Genome()

    def barf_session(self, line):
        if len(line) == 0:
            print("Usage: barfsession <directory>\n")
            return

        os.system('mkdir '+line)
        
        # Write the gff
        with open(line+'/gag.gff', 'w') as gff:
            for gene in self.genome.gff.genes:
                gff.write(gene.to_gff())

        # Write the fasta
        with open(line+'/gag.fasta', 'w') as fasta:
            fasta.write(self.genome.fasta.write_string())

        # Write the annotations
        self.genome.annot.write_to_file(line+'/gag.trinotate')

    def load_session(self, line):
        # Read the gff
        self.read_gff(line+'/gag.gff')

        # Read the fasta
        self.read_fasta(line+'/gag.fasta')

        # Read the annotations
        self.read_trinotate(line+'/gag.trinotate')


## Reading in files

    def read_fasta(self, line):
        self.genome.fasta = Fasta()
        self.genome.fasta.read_file(line)

    def read_gff(self, line):
        self.genome.gff = GFF()
        with open(line, 'rb') as gfffile:
            gffreader = csv.reader(gfffile, delimiter='\t')
            self.genome.gff.read_file(gffreader)
        return self.genome.gff

    def read_trinotate(self, line):
        self.genome.annot = Annotator()
        self.genome.annot.read_from_file(line)


## Manipulate genome

    def apply_bed(self, line):
        bed = Bed()
        with open(line, 'rb') as bedfile:
            bedreader = csv.reader(bedfile, delimiter='\t')
            bed.read_file(bedreader)
            self.genome.fasta.apply_bed(bed)
            self.genome.gff.apply_bed(bed)
            self.genome.gff.remove_empty_genes()

    def duct_tape_seq_frames(self, line):
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


## Output info to console

    def barf_gff(self, line):
        for gene in self.genome.gff.genes:
            if gene.name == line:
                print(gene.to_gff())

    def barf_seq(self, line):
        args = line.split(' ')
        print(str(self.genome.fasta.get_subseq(args[0], [int(args[1]), int(args[2])]))+'\n')

    def barf_gene_tbl(self, line):
        self.genome.write_file(sys.stdout, set(line.split()))

    def barf_err_subset(self, line):
        args = line.split()
        if len(args) < 2:
            print("Usage: barferrsubset <directory> <errorcode>\n")
            return
        outdir = args[0]
        err = args[1]
        # TODO ?

## Output info to file

    def write_tbl(self, line):
        with open(line, 'w') as outFile:
            self.genome.write_file(outFile)






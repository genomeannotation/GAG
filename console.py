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
from src.feature_tbl import FeatureTbl
from src.annotator import Annotator

class GagCmd(cmd.Cmd):

    def __init__(self):
        cmd.Cmd.__init__(self)
        self.prompt = "GAG> "

        readline.set_history_length(1000)
        readline.read_history_file('.gaghistory')

        self.feature_tbl = FeatureTbl()
                

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_readgff(self):
        print("Usage: readgff <file_name>\n")
        print("Read the gff file. Any unsaved changes")
        print("to the currently loaded gff will be lost.\n")

    def do_readgff(self, line):
        self.feature_tbl.gff = GFF()
        with open(line, 'rb') as gfffile:
            gffreader = csv.reader(gfffile, delimiter='\t')
            self.feature_tbl.gff.read_file(gffreader)

    def help_readfasta(self):
        print("Usage: readfasta <file_name>\n")
        print("Read the fasta file. Any unsaved changes")
        print("to the currently loaded fasta will be lost.\n")

    def do_readfasta(self, line):
        self.feature_tbl.fasta = Fasta()
        self.feature_tbl.fasta.read_file(line)

    def help_readtrinotate(self):
        print("Usage: readtrinotate <file_name>\n")
        print("Read annotations from a trinotate file.\n")

    def do_readtrinotate(self, line):
        self.feature_tbl.annot = Annotator()
        self.feature_tbl.annot.read_from_file(line)

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
            self.feature_tbl.fasta.apply_bed(bed)
            self.feature_tbl.gff.apply_bed(bed)
            self.feature_tbl.gff.remove_empty_genes()

    def help_writetbl(self):
        print("Usage: writetbl <file_name>\n")
        print("Write a sweet feature table to the specified file.\n")

    def do_writetbl(self, line):
        with open(line, 'w') as outFile:
            self.feature_tbl.write_file(outFile)
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
            for gene in self.feature_tbl.gff.genes:
                gff.write(gene.to_gff())
            gff.close()

        # Write the fasta
        with open(line+'/gag.fasta', 'w') as fasta:
            fasta.write(self.feature_tbl.fasta.write_string())
            fasta.close()

        # Write the annotations
        self.feature_tbl.annot.write_to_file(line+'/gag.trinotate')

    def do_loadsession(self, line):
        # Read the gff
        self.do_readgff(line+'/gag.gff')

        # Read the fasta
        self.do_readfasta(line+'/gag.fasta')

        # Read the annotations
        self.do_readtrinotate(line+'/gag.trinotate')

    def do_barfgenegff(self, line):
        for gene in self.feature_tbl.gff.genes:
            if gene.name == line:
                print(gene.to_gff())

    def do_barfgenetbl(self, line):
        self.feature_tbl.write_file(sys.stdout, set(line.split()))

    def help_exit(self):
        print("Exit this console.\n")

    def do_exit(self, line):
        return True

if __name__ == '__main__':
    GagCmd().cmdloop("Welcome to the GAG console!")

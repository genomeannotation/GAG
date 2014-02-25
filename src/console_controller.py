#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import csv
import subprocess
import glob
from src.fasta import Fasta
from src.gff import GFF
from src.genome import Genome
from src.annotator import Annotator
from src.translate import translate

def foo():
    return 3

class ConsoleController:

## Setup, loading and saving sessions

    def __init__(self):
        self.genome = Genome()
        self.input = ''

    def barf_folder(self, line):
        if len(line) == 0:
            print("Usage: barffolder <directory>\n")
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
        
    def load_folder(self, line):
        # Get filenames
        gffs = glob.glob(line + '/*.gff')
        fastas = glob.glob(line + '/*.fasta')
        trinotates = glob.glob(line + '/*.trinotate')

        # Read the gff
        if gffs:
            self.read_gff(gffs[0])
        else:
            sys.stderr.write("Couldn't find .gff file in " + line + "\n")
            return

        # Read the fasta
        if fastas:
            self.read_fasta(fastas[0])
        else:
            sys.stderr.write("Couldn't find .fasta file in " + line + "\n")
            return

        # Read the annotations
        if trinotates:
            self.read_trinotate(line+'/gag.trinotate')
        else:
            sys.stderr.write("Did not find .trinotate file; no functional annotations available.\n")


    def ls(self, line):
        proc = subprocess.Popen(['ls '+line], stdout=subprocess.PIPE, \
                stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

    def cat(self, line):
        proc = subprocess.Popen(['cat '+line], stdout=subprocess.PIPE, \
                stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

    def grep(self, line):
        proc = subprocess.Popen(['grep '+line], stdout=subprocess.PIPE, \
                stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

    def sed(self, line):
        proc = subprocess.Popen(['sed '+line], stdout=subprocess.PIPE, \
                stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

    def sort(self, line):
        proc = subprocess.Popen(['sort '+line], stdout=subprocess.PIPE, \
                stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

    def uniq(self, line):
        proc = subprocess.Popen(['uniq '+line], stdout=subprocess.PIPE, \
                stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

    def barf(self, line):
        proc = subprocess.Popen(['echo '+line], stdout=subprocess.PIPE, \
                stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

## Assorted utilities

    def status(self):
        result = ""
        if self.genome.fasta:
            result += "Fasta: " + str(self.genome.fasta)
        else:
            result += "Fasta: no fasta\n"
        if self.genome.gff:
            result += "GFF: " + str(self.genome.gff)
        else:
            result += "GFF: no gff\n"
        return result

    def barftofile(self, line):
        args = line.split()

        with open(args[0], 'w') as f:
            if len(args) > 1:
                for arg in args[1:]:
                    f.write(arg+' ')
            else:
                f.write(self.input)


## Reading in files

    def read_fasta(self, line):
        self.genome.fasta = Fasta()
        self.genome.fasta.read(open(line, 'r'))

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

    def ducttape(self):
        min_first_cds_segment_length = 3
        min_cds_length = 150
        if self.genome.gff:
            self.genome.rename_maker_mrnas()
            self.genome.gff.remove_first_cds_segment_if_shorter_than(min_first_cds_segment_length)
            self.genome.create_starts_and_stops() 
            self.genome.gff.remove_mrnas_with_cds_shorter_than(min_cds_length)

    def create_starts_and_stops(self):
        self.genome.create_starts_and_stops() 

    def subset_genome(self, line):
        # parse args
        args = line.split()
        if args:
            self.genome.fasta.subset_fasta(args)
            self.genome.gff.subset_gff(args)

    def subset_fasta(self):
        # line parameter is not used, but Cmd likes to pass it so there it is.
        self.genome.fasta.subset_fasta(self.seqlist)

    def subset_gff(self):
        self.genome.gff.subset_gff(self.seqlist)

    def duct_tape_seq_frames(self, line):
        result = ''
        args = None        

        if len(line) > 0:
            args = line.split()
        else:
            args = self.input.split('\n')

        for yarg in args: # I'm a pirate
            result += self.ducttape_mrna_seq_frame(yarg)
        return result


    def removemrna(self, line):
        args = None        

        if len(line) > 0:
            args = line.split()
        else:
            args = self.input.split('\n')

        for name in args:
            eraseGenes = []
            for gene in self.genome.gff.genes:
                erase = []
                for mrna in gene.mrnas:
                    if mrna.name == name:
                        erase.append(mrna)
                for mrna in erase:
                    gene.mrnas.remove(mrna)
                if len(gene.mrnas) == 0:
                    eraseGenes.append(gene)
            for gene in eraseGenes:
                self.genome.gff.genes.remove(gene)
                        
    def remove_genes_marked_for_removal(self, line):
        self.genome.remove_genes_marked_for_removal()

    def rename_maker_mrnas(self):
        self.genome.rename_maker_mrnas()

    def ducttape_mrna_seq_frame(self, name):
        for gene in self.genome.gff.genes:
            for mrna in gene.mrnas:
                if mrna.name == name:
                    seq = self.genome.fasta.get_subseq(gene.seq_name, \
                            [mrna.cds.indices[0]]) #first segment
                    if seq == None:
                        return "Failed to fix "+name+\
                               ": sequence does not exist.\n" 
                    elif len(seq) < 6:
                        return "Failed to fix "+name+\
                               ": sequence less than 6 base pairs.\n"

                    pseq1 = translate(seq, 1, '+')
                    pseq2 = translate(seq, 2, '+')
                    pseq3 = translate(seq, 3, '+')
                    nseq1 = translate(seq, 1, '-')
                    nseq2 = translate(seq, 2, '-')
                    nseq3 = translate(seq, 3, '-')

                    annotEntry = self.genome.annot.get_entry(name)
                    if annotEntry:
                        pepSeq = annotEntry[9]
                        if pepSeq == None:
                            return "Failed to fix "+name+\
                                   ": trinotate missing peptide sequence.\n"

                        oldphase = mrna.cds.phase[0]
                        if pseq1 and pepSeq.find(pseq1[:-1]) == 0:
                            gene.strand = '+'
                            mrna.cds.phase[0] = 0
                        elif pseq2 and pepSeq.find(pseq2[:-1]) == 0:
                            gene.strand = '+'
                            mrna.cds.phase[0] = 1
                        elif pseq3 and pepSeq.find(pseq3[:-1]) == 0:
                            gene.strand = '+'
                            mrna.cds.phase[0] = 2
                        elif nseq1 and pepSeq.find(nseq1[:-1]) == 0:
                            gene.strand = '-'
                            mrna.cds.phase[0] = 0
                        elif nseq2 and pepSeq.find(nseq2[:-1]) == 0:
                            gene.strand = '-'
                            mrna.cds.phase[0] = 1
                        elif nseq3 and pepSeq.find(nseq3[:-1]) == 0:
                            gene.strand = '-'
                            mrna.cds.phase[0] = 2
                        else:
                            return "Failed to fix "+name+\
                                   ": no matching translation.\n"
                        return "Fixed "+name+" from phase "+str(oldphase)+\
                               " to phase "+str(mrna.cds.phase[0])+"\n"
                    else:
                        return "Failed to fix "+name+\
                               ": trinotate entry doesn't exist.\n"
        return "Failed to fix "+name+": mRNA doesn't exist.\n"

    def remove_gene(self, line):
        args = []
        if len(line) > 0:
            args = line.split()
        else:
            args = self.input.split('\n')
        for arg in args:
            self.genome.remove_all_gene_segments(arg)

    def remove_mrnas_with_cds_shorter_than(self, line):
        min_length = int(line)
        self.genome.remove_mrnas_with_cds_shorter_than(min_length)

    def trim_region(self, line):
        args = []
        if len(line) > 0:
            args = line.split()
            if len(args) != 3:
                sys.stderr.write("Error: ConsoleController.trim_region \
                                  requires 3 args\n")
            else:
                seq = args[0]
                start = int(args[1])
                stop = int(args[2])
                self.genome.trim_region(seq, start, stop)
        else:
            lines = self.input.split('\n')
            for entry in lines:
                entries = entry.split()
                if len(entries) != 3:
                    sys.stderr.write("Error: ConsoleController.trim_region " +
                                      "requires 3 args\n")
                    sys.stderr.write("This was the input: " + entry + "\n")
                    sys.stderr.write("Moving on to next input...\n")
                else:
                    seq = entries[0]
                    start = int(entries[1])
                    stop = int(entries[2])
                    self.genome.trim_region(seq, start, stop)

    def remove_seq(self, line):
        args = []
        if len(line) > 0:
            args = line.split()
            seq_id = args[0]
            if len(args) == 2 and args[1] == '-F':
                self.genome.remove_seq(seq_id, force=True)
            else:
                self.genome.remove_seq(seq_id)

    def check_gene_for_invalid_begin_or_end(self, line):
        args = []
        if len(line) > 0:
            args = line.split()
        else:
            args = self.input.split('\n')
        for arg in args:
            self.genome.check_gene_for_invalid_begin_or_end(arg)

    def invalidate_region(self, line):
        if len(line) > 0:
            args = line.split()
            self.genome.invalidate_region(args[0], int(args[1]), int(args[2]))
        else:
            lines = self.input.split('\n')
            for line in lines:
                args = line.split()
                if not args:
                    continue
                self.genome.invalidate_region(args[0], \
                        int(args[1]), int(args[2]))


## Output info to console

    def barf_gff(self, line):
        for gene in self.genome.gff.genes:
            if gene.name == line:
                return gene.to_gff()

    def barf_seq(self, line):
        args = line.split(' ')
        return str(self.genome.fasta.get_subseq(args[0], \
                [[int(args[1]), int(args[2])]]))+'\n'

    def barf_cds_seq(self, line):
        name = line

        for gene in self.genome.gff.genes:
            for mrna in gene.mrnas:
                if mrna.name == name and mrna.cds:
                    return mrna.cds.extract_sequence(self.genome.fasta, \
                            gene.seq_name, gene.strand)

        return "Error: Couldn't find mRNA.\n"

    def barf_gene_tbl(self, line):
        return self.genome.write_string(set(line.split()))

## Output info to file

    def write_tbl(self, line):
        with open(line, 'w') as outFile:
            outFile.write(self.genome.write_string())
            outFile.close()

    def write_fasta(self, line):
        with open(line, 'w') as outFile:
            outFile.write(self.genome.fasta.write_string())


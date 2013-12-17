#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import csv
import subprocess
import StringIO
import pdb
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
        self.input = ''
        self.tbl2asn_executable = None
        self.template_file = None
        self.seqlist = []

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

    def ls(self, line):
        proc = subprocess.Popen(['ls '+line], stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

    def cat(self, line):
        proc = subprocess.Popen(['cat '+line], stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

    def grep(self, line):
        proc = subprocess.Popen(['grep '+line], stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

    def sed(self, line):
        proc = subprocess.Popen(['sed '+line], stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

    def sort(self, line):
        proc = subprocess.Popen(['sort '+line], stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

    def uniq(self, line):
        proc = subprocess.Popen(['uniq '+line], stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

    def barf(self, line):
        proc = subprocess.Popen(['echo '+line], stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

## Assorted utilities

    def add_seq(self, line):
        args = []
        if len(line) > 0:
            args = line.split()
        else:
            args = self.input.split('\n')

        for arg in args:
            self.seqlist.append(arg)

    def clear_seqlist(self):
        del(self.seqlist[:])

    def add_template_file(self, line):
        self.template_file = line

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
        if self.template_file:
            result += "Template File: " + self.template_file + "\n"
        else:
            result += "Template File: no template file\n"
        if self.seqlist:
            result += "Seqlist: " + str(self.seqlist) + "\n"
        else:
            result += "Seqlist: no seqlist\n"
        if self.tbl2asn_executable:
            result += "Tbl2asn Executable: " + self.tbl2asn_executable + "\n"
        else:
            result += "Tbl2asn Executable: no tbl2asn executable\n"
        return result

    def barftofile(self, line):
        args = line.split()

        with open(args[0], 'w') as f:
            if len(args) > 1:
                for arg in args[1:]:
                    f.write(arg+' ')
            else:
                f.write(self.input)

    def barffromfile(self, line):
        result = ''
        with open(line.strip(), 'r') as f:
            for fline in f:
                result += fline
        return result
            


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
            if self.ducttape_mrna_seq_frame(yarg):
                result += 'Fixed '+yarg+'.\n'
            else:
                result += 'Failed to fix '+yarg+'.\n'
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
                        

    def ducttape_mrna_seq_frame(self, name):
        for gene in self.genome.gff.genes:
            for mrna in gene.mrnas:
                if mrna.name == name:
                    seq = self.genome.fasta.get_subseq(gene.seq_name, mrna.cds.indices[0])
                    if seq == None or len(seq) < 6:
                        return False
                    
                    pseq1 = translate(seq, 1, '+')[:-1]
                    pseq2 = translate(seq, 2, '+')[:-1]
                    pseq3 = translate(seq, 3, '+')[:-1]
                    nseq1 = translate(seq, 1, '-')[:-1]
                    nseq2 = translate(seq, 2, '-')[:-1]
                    nseq3 = translate(seq, 3, '-')[:-1]

                    annotEntry = self.genome.annot.get_entry(name)
                    if annotEntry:
                        pepSeq = annotEntry[9]
                        if pepSeq == None:
                            return False

                        if pseq1 and pepSeq.find(pseq1) != -1:
                            gene.strand = '+'
                            mrna.cds.phase[0] = 0
                        elif pseq2 and pepSeq.find(pseq2) != -1:
                            gene.strand = '+'
                            mrna.cds.phase[0] = 1
                        elif pseq3 and pepSeq.find(pseq3) != -1:
                            gene.strand = '+'
                            mrna.cds.phase[0] = 2
                        elif nseq1 and pepSeq.find(nseq1) != -1:
                            gene.strand = '-'
                            mrna.cds.phase[0] = 0
                        elif nseq2 and pepSeq.find(nseq2) != -1:
                            gene.strand = '-'
                            mrna.cds.phase[0] = 1
                        elif nseq3 and pepSeq.find(nseq3) != -1:
                            gene.strand = '-'
                            mrna.cds.phase[0] = 2
                        else:
                            return False
                        return True
                    else:
                        return False
        return False

    def remove_all_gene_segments(self, line):
        args = []
        if len(line) > 0:
            args = line.split()
        else:
            args = self.input.split('\n')
        for arg in args:
            self.genome.remove_all_gene_segments(arg)
        

## Output info to console

    def barf_gff(self, line):
        for gene in self.genome.gff.genes:
            if gene.name == line:
                return gene.to_gff()

    def barf_seq(self, line):
        args = line.split(' ')
        return str(self.genome.fasta.get_subseq(args[0], [int(args[1]), int(args[2])]))+'\n'

    def barf_gene_tbl(self, line):
        return self.genome.write_string(set(line.split()))

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
            outFile.write(self.genome.write_string())
            outFile.close()

    def write_fasta(self, line):
        with open(line, 'w') as outFile:
            outFile.write(self.genome.fasta.write_string())

## New, Exciting tbl2asn integration
    def set_tbl2asn_executable(self, line):
        self.tbl2asn_executable = line

    def prep_tbl2asn(self, line):
        if os.path.exists(line):
            sys.stderr.write("Sorry, looks like " + line + " already exists.\n")
            sys.stderr.write("Please try command again with another directory name.\n")
            return
        elif not self.template_file:
            sys.stderr.write("No template file specified. Try 'addtemplatefile'\n")
            return
        else:
            # create tbl2asn directory
            os.system('mkdir ' + line)
            # symlink template file 
            template_abs_path = os.path.abspath(self.template_file)
            os.system('ln -s ' + template_abs_path + ' ' + line + '/gag.sbt')
            # write fasta file
            self.write_fasta(line + '/gag.fsa')
            # write tbl file
            self.write_tbl(line + "/gag.tbl")

    def ready_for_tbl2asn(self, line):
        if not self.tbl2asn_executable:
            return False
        elif not os.path.isdir(line):
            return False
        elif not os.path.exists(line + "/gag.fsa"):
            return False
        elif not os.path.exists(line + "/gag.tbl"):
            return False
        elif not os.path.exists(line + "/gag.sbt"):
            return False
        else:
            return True

    def run_tbl2asn(self, line):
        if self.ready_for_tbl2asn(line):
            tbl2asn_command = self.tbl2asn_executable + " -p " + line
            tbl2asn_command += ' -j "[organism=Bactrocera dorsalis][tech=WGS]" -M b -V vb -c f -Z discrep'
            os.system(tbl2asn_command)
        else:
            sys.stderr.write("Sorry, unable to run tbl2asn in " + line + ". Try prep_tbl2asn or settbl2asnexecutable first.")
            






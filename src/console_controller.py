#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import csv
import subprocess
import glob
import copy
from src.fasta_reader import FastaReader
from src.gff_reader import GFFReader
from src.annotator import Annotator
from src.filter_manager import FilterManager
from src.stats_manager import StatsManager
from src.translate import translate

class ConsoleController:

## Setup, loading and saving sessions

    def __init__(self):
        self.seqs = []
        self.annot = Annotator()
        self.filter_mgr = FilterManager()
        self.stats_mgr = StatsManager()
        self.input = ''

    def barf_folder(self, line):
        if len(line) == 0:
            sys.stderr.write("Usage: barffolder <directory>\n")
            return

        os.system('mkdir '+line)
        
        # Write the gff
        with open(line+'/gag.gff', 'w') as gff:
            for seq in self.seqs:
                gff.write(seq.to_gff())

        # Write the fasta
        with open(line+'/gag.fasta', 'w') as fasta:
            for seq in self.seqs:
                fasta.write(seq.to_fasta())

        # Write the annotations
        self.annot.write_to_file(line+'/gag.trinotate')
        
    def load_folder(self, line):
        if not line:
            line = "."
        # Get filenames
        sys.stderr.write("Locating files...\n")
        gffs = glob.glob(line + '/*.gff')
        fastas = glob.glob(line + '/*.fasta')
        trinotates = glob.glob(line + '/*.trinotate')

        # Make sure there's only one of each file type
        if len(gffs) > 1:
            sys.stderr.write("Found more than one gff file; no genome loaded.\n")
            return
        elif len(fastas) > 1:
            sys.stderr.write("Found more than one fasta file; no genome loaded.\n")
            return
        elif len(trinotates) > 1:
            sys.stderr.write("Found more than one trinotate file; no genome loaded.\n")
            return

        # Read the fasta
        if fastas:
            sys.stderr.write("Reading fasta...\n")
            self.read_fasta(fastas[0])
            sys.stderr.write("Done.\n")
        else:
            sys.stderr.write("Couldn't find .fasta file in " + line + "\n")
            return

        # Read the gff
        if gffs:
            sys.stderr.write("Reading gff...\n")
            self.read_gff(gffs[0])
            sys.stderr.write("Done.\n")
        else:
            sys.stderr.write("Couldn't find .gff file in " + line + "\n")
            return

        # Read the annotations
        if trinotates:
            sys.stderr.write("Reading trinotate...\n")
            self.read_trinotate(line+'/gag.trinotate')
            sys.stderr.write("Done.\n")
        else:
            sys.stderr.write("Did not find .trinotate file; no functional annotations available.\n")

        # Read in stats
        sys.stderr.write("Calculating statistics on genome...\n")
        for seq in self.seqs:
            self.stats_mgr.update_ref(seq.stats())
        sys.stderr.write("Done.\n\n")

    def set_filter_arg(self, filter_name, filter_arg, val):
        self.filter_mgr.set_filter_arg(filter_name, filter_arg, val)
        
    def get_filter_arg(self, filter_name, filter_arg):
        return str(self.filter_mgr.get_filter_arg(filter_name, filter_arg))
        
    def get_filter_help(self, filter_name):
        return self.filter_mgr.get_filter_help(filter_name)
        
    def apply_filters(self):
        for seq in self.seqs:
            self.filter_mgr.apply_filters(seq)

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
        return "Number of seqs: " + str(len(self.seqs))

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
        reader = FastaReader()
        self.seqs = reader.read(open(line, 'r'))

    def read_gff(self, line):
        gffreader = GFFReader()
        reader = open(line, 'rb')
        genes = gffreader.read_file(reader)
        for gene in genes:
            self.add_gene(gene)

    def read_trinotate(self, line):
        self.annot.read_from_file(line)



## Manipulate genome

    def ducttape(self):
        # TODO
        pass
        #min_first_cds_segment_length = 3
        #min_cds_length = 150
        #if self.genome.genes:
            #self.genome.rename_maker_mrnas()
            #self.genome.remove_first_cds_segment_if_shorter_than(min_first_cds_segment_length)
            #self.genome.create_starts_and_stops() 
            #self.genome.remove_mrnas_with_cds_shorter_than(min_cds_length)

    def create_starts_and_stops(self):
        #TODO pleeeeeease pass seq.bases instead of whole seq.
        for seq in self.seqs:
            seq.create_starts_and_stops()

    def subset_genome(self, line):
        # parse args
        args = line.split()
        if args:
            self.seqs = [s for s in self.seqs if s.header in args]

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

        for seq in self.seqs:
            seq.remove_mrna(args)
                        
    def rename_maker_mrnas(self):
        locus_tag = self.get_locus_tag()
        count = 1000000
        for seq in self.seqs:
            for gene in seq.genes:
                for mrna in gene.mrnas:
                    if mrna.is_maker_mrna():
                        old_name = mrna.identifier
                        new_name = locus_tag + '_' + str(count)
                        mrna.identifier = new_name
                        self.annot.rename_mrna(old_name, new_name)
                        count += 1

    def ducttape_mrna_seq_frame(self, name):
        for seq in self.seqs:
            for gene in seq.genes:
                for mrna in gene.mrnas:
                    if mrna.identifier == name:
                        subseq = seq.get_subseq(mrna.cds.indices[0], mrna.cds.indices[1]) #first segment
                        if subseq == None:
                            return "Failed to fix "+name+\
                                   ": sequence does not exist.\n" 
                        elif len(subseq) < 6:
                            return "Failed to fix "+name+\
                                   ": sequence less than 6 base pairs.\n"

                        pseq1 = translate(subseq, 1, '+')
                        pseq2 = translate(subseq, 2, '+')
                        pseq3 = translate(subseq, 3, '+')
                        nseq1 = translate(subseq, 1, '-')
                        nseq2 = translate(subseq, 2, '-')
                        nseq3 = translate(subseq, 3, '-')

                        annotEntry = self.annot.get_entry(name)
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
        for seq in self.seqs:
            seq.remove_gene(args)

    def remove_mrnas_with_cds_shorter_than(self, line):
        min_length = int(line)
        for seq in self.seqs:
            seq.remove_mrnas_with_cds_shorter_than(min_length)

    def trim_region(self, line):
        args = []
        if len(line) > 0:
            args = line.split()
            if len(args) != 3:
                sys.stderr.write("Error: ConsoleController.trim_region \
                                  requires 3 args\n")
            else:
                seq_name = args[0]
                start = int(args[1])
                stop = int(args[2])
                for seq in self.seqs:
                    if seq.header == seq_name:
                        seq.trim_region(start, stop)
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
                    # TODO too many loops, could be nicer
                    seq_name = entries[0]
                    start = int(entries[1])
                    stop = int(entries[2])
                    for seq in self.seqs:
                        if seq.header == seq_name:
                            seq.trim_region(start, stop)

    def remove_seq(self, line):
        # TODO take multiple args?
        self.seqs = [s for s in self.seqs if s.header != line]

    def check_gene_for_invalid_begin_or_end(self, line):
        # TODO this should probably just check all genes instead of taking args
        #args = []
        #if len(line) > 0:
            #args = line.split()
        #else:
            #args = self.input.split('\n')
        #for arg in args:
            #self.genome.check_gene_for_invalid_begin_or_end(arg)
        pass

    def invalidate_region(self, line):
        # TODO return error messages on invalid args
        # TODO not working
        if len(line) > 0:
            args = line.split()
            seq_name = args[0]
            start = args[1]
            stop = args[2]
            for seq in self.seqs:
                if seq.header == seq_name:
                    seq.invalidate_region(start, stop)
        else:
            lines = self.input.split('\n')
            for line in lines:
                args = line.split()
                if not args:
                    continue
                seq_name = args[0]
                start = args[1]
                stop = args[2]
                for seq in self.seqs:
                    if seq.header == seq_name:
                        seq.invalidate_region(start, stop)


## Output info to console

    def barf_gene_gff(self, line):
        for seq in self.seqs:
            if seq.contains_gene(line):
                return seq.gene_to_gff(line)

    def barf_seq(self, line):
        args = line.split(' ')
        if len(args) != 3:
            return "Usage: barfseq <seq_id> <start_index> <end_index>\n"
        seq_id = args[0]
        start = int(args[1])
        stop = int(args[2])
        for seq in self.seqs:
            if seq.header == seq_id:
                return seq.get_subseq(start, stop)

    def barf_cds_seq(self, line):
        name = line
        for seq in self.seqs:
            if seq.contains_mrna(name):
                return seq.extract_cds_seq(name)
        return "Error: Couldn't find mRNA.\n"

    def barf_gene_tbl(self, line):
        # TODO this used to take multiple gene_ids? but do we care?
        output = ">Feature SeqId\n"
        for seq in self.seqs:
            if seq.contains_gene(line):
                output += seq.gene_to_tbl(line)
        return output

    def stats(self):
        self.stats_mgr.clear_alt()
        sys.stderr.write("Calculating statistics on genome...\n")
        for seq in self.seqs:
            cseq = copy.deepcopy(seq)
            self.filter_mgr.apply_filters(cseq)
            self.stats_mgr.update_alt(cseq.stats())
        return self.stats_mgr.summary()

## Output info to file

    def write_tbl(self, line):
        if os.path.exists(line):
            return line + " already exists; please try another filename\n"
        with open(line, 'w') as outFile:
            outFile.write(">Feature SeqId\n")
            for seq in self.seqs:
                outFile.write(seq.to_tbl())
            outFile.close()
        return ".tbl file written to " + line + "\n"

    def write_fasta(self, line):
        with open(line, 'w') as outFile:
            for seq in self.seqs:
                outFile.write(seq.to_fasta())

## Utilities

    def add_gene(self, gene):
        for seq in self.seqs:
            if seq.header == gene.seq_name:
                seq.add_gene(gene)

    def get_locus_tag(self):
        locus_tag = ""
        for seq in self.seqs:
            if locus_tag:
                break
            else:
                locus_tag = seq.get_locus_tag()
        return locus_tag
    

#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import subprocess
import glob
import copy
from src.fasta_reader import FastaReader
from src.gff_reader import GFFReader
from src.annotator import Annotator
from src.filter_manager import FilterManager
from src.stats_manager import StatsManager

class ConsoleController:

    no_genome_message = "It looks like no genome is currently loaded. Try the 'loadfolder' command.\n"+\
            "Type 'help loadfolder' to learn how to use it, or just 'help' for general advice.\n"

## Setup, loading and saving sessions

    def __init__(self):
        self.seqs = []
        self.annot = Annotator()
        self.filter_mgr = FilterManager()
        self.stats_mgr = StatsManager()
        self.input = ''

    def barf_folder(self, line):
        if not self.seqs:
            return self.no_genome_message
        elif len(line) == 0:
            sys.stderr.write("Usage: barffolder <directory>\n")
            return
        else:
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
            return "Genome written to " + line
        
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

        # Clear stats; read in new stats; display stats on reference genome
        self.stats_mgr.clear_all()
        for seq in self.seqs:
            self.stats_mgr.update_ref(seq.stats())
        sys.stderr.write("Done.\n\n")
        print(self.stats())

    def set_filter_arg(self, filter_name, filter_arg, val):
        self.filter_mgr.set_filter_arg(filter_name, filter_arg, val)
        
    def get_filter_arg(self, filter_name, filter_arg):
        return str(self.filter_mgr.get_filter_arg(filter_name, filter_arg))
        
    def get_filter_help(self, filter_name):
        return self.filter_mgr.get_filter_help(filter_name)
        
    def apply_filters(self):
        for seq in self.seqs:
            self.filter_mgr.apply_filters(seq)

    def barf(self, line):
        proc = subprocess.Popen(['echo '+line], stdout=subprocess.PIPE, \
                stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(self.input)
        return out

## Assorted utilities

    def status(self):
        if not self.seqs:
            return self.no_genome_message
        else:
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

    def create_starts_and_stops(self):
        if not self.seqs:
            return self.no_genome_message
        else:
            #TODO pleeeeeease pass seq.bases instead of whole seq.
            for seq in self.seqs:
                seq.create_starts_and_stops()

    def subset_genome(self, line):
        if not self.seqs:
            return self.no_genome_message
        else:
            # parse args
            args = line.split()
            if args:
                self.seqs = [s for s in self.seqs if s.header in args]

    def removemrna(self, line):
        args = None        

        if len(line) > 0:
            args = line.split()
        else:
            args = self.input.split('\n')

        for seq in self.seqs:
            seq.remove_mrna(args)
                        
    def remove_gene(self, line):
        if not self.seqs:
            return self.no_genome_message
        else:
            args = []
            if len(line) > 0:
                args = line.split()
            else:
                args = self.input.split('\n')
            for seq in self.seqs:
                seq.remove_gene(args)

    def trim_region(self, line):
        # TODO this is hideous 
        if not self.seqs:
            return self.no_genome_message
        else:
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
        if not self.seqs:
            return self.no_genome_message
        else:
            # TODO take multiple args?
            self.seqs = [s for s in self.seqs if s.header != line]

    def invalidate_region(self, line):
        if not self.seqs:
            return self.no_genome_message
        else:
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
        if not self.seqs:
            return self.no_genome_message
        else:
            for seq in self.seqs:
                if seq.contains_gene(line):
                    return seq.gene_to_gff(line)

    def barf_seq(self, line):
        if not self.seqs:
            return self.no_genome_message
        else:
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
        if not self.seqs:
            return self.no_genome_message
        else:
            name = line
            for seq in self.seqs:
                if seq.contains_mrna(name):
                    return seq.extract_cds_seq(name)
            return "Error: Couldn't find mRNA.\n"

    def barf_gene_tbl(self, line):
        if not self.seqs:
            return self.no_genome_message
        else:
            # TODO this used to take multiple gene_ids? but do we care?
            output = ">Feature SeqId\n"
            for seq in self.seqs:
                if seq.contains_gene(line):
                    output += seq.gene_to_tbl(line)
            return output

    def stats(self):
        if not self.seqs:
            return self.no_genome_message
        else:
            first_line = "Number of sequences:   " + str(len(self.seqs)) + "\n"
            if self.filter_mgr.dirty:
                self.stats_mgr.clear_alt()
                sys.stderr.write("Calculating statistics on genome...\n")
                for seq in self.seqs:
                    cseq = copy.deepcopy(seq)
                    self.filter_mgr.apply_filters(cseq)
                    print("updating alt")
                    self.stats_mgr.update_alt(cseq.stats())
                self.filter_mgr.dirty = False
            return first_line + self.stats_mgr.summary()

## Output info to file

    def write_tbl(self, line):
        if not self.seqs:
            return self.no_genome_message
        elif not line:
            return "Usage: writetbl <filename>\n"
        else:
            if os.path.exists(line):
                return line + " already exists; please try another filename\n"
            with open(line, 'w') as out_file:
                out_file.write(">Feature SeqId\n")
                for seq in self.seqs:
                    out_file.write(seq.to_tbl())
                out_file.close()
            return ".tbl file written to " + line + "\n"

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
    
    def clear_seqs(self):
        self.seqs[:] = []

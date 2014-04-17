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
from src.seq_fixer import SeqFixer

class ConsoleController:

    no_genome_message = "It looks like no genome is currently loaded. Try the 'load' command.\n"+\
            "Type 'help load' to learn how to use it, or just 'help' for general advice.\n"

## Setup, loading and saving sessions

    def __init__(self):
        self.seqs = []
        self.annot = Annotator()
        self.filter_mgr = FilterManager()
        self.stats_mgr = StatsManager()
        self.seq_fixer = SeqFixer()

    def genome_is_loaded(self):
        if self.seqs:
            return True
        else:
            return False

    def barf_folder(self, line):
        if not self.seqs:
            return self.no_genome_message
        elif len(line) == 0:
            sys.stderr.write("Usage: barffolder <directory>\n")
            return
        else:
            # Create directory, open files
            os.system('mkdir '+line)
            gff = open(line+'/gag.gff', 'w')
            tbl = open(line+'/gag.tbl', 'w')
            fasta = open(line+'/gag.fasta', 'w')

            # Deep copy each seq, apply fixes and filters, write
            sys.stderr.write("Writing gff, tbl and fasta...\n")
            for seq in self.seqs:
                cseq = copy.deepcopy(seq)
                self.seq_fixer.fix(seq)
                self.filter_mgr.apply_filters(cseq)
                gff.write(seq.to_gff())
                tbl.write(seq.to_tbl())
                fasta.write(seq.to_fasta())

            # Close files
            gff.close()
            tbl.close()
            fasta.close()

            # Write the annotations
            sys.stderr.write("Writing trinotate...\n")
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
        sys.stderr.write("Genome loaded.\n\n")
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

    def fix_terminal_ns(self):
        self.seq_fixer.fix_terminal_ns()
        return "Terminal Ns will now be fixed."

    def fix_start_stop_codons(self):
        self.seq_fixer.fix_start_stop_codons()
        return "Will verify and create start/stop codons."

## Assorted utilities

    def get_n_seq_ids(self, number):
        """Returns a message indicating the first n seq_ids in the genome.

        If no seqs loaded, returns a message to that effect. If fewer than n
        seqs loaded, returns the seq_ids of those seqs."""
        if not self.seqs:
            return "No sequences currently in memory.\n"
        else:
            if len(self.seqs) < number:
                number = len(self.seqs)
            seq_list = []
            for seq in self.seqs:
                seq_list.append(seq.header)
                if len(seq_list) == number:
                    break
            result = "First " + str(len(seq_list)) + " seq ids are: "
            result += format_list_with_strings(seq_list)
            return result

    def get_n_gene_ids(self, number):
        """Returns a message indicating the first n gene_ids in the genome.

        If no genes are present, returns a message to that effect. If fewer than n
        genes are loaded, returns the gene_ids of those genes."""
        genes_list = []
        while len(genes_list) < number:
            for seq in self.seqs:
                genes_list.extend(seq.get_gene_ids())
        # List may now contain more than 'number' ids, or it may contain zero
        if not genes_list:
            return "No genes currently in memory.\n"
        if len(genes_list) > number:
            genes_list = genes_list[:number]
        result = "First " + str(len(genes_list)) + " gene ids are: "
        result += format_list_with_strings(genes_list)
        return result

    def get_n_mrna_ids(self, number):
        """Returns a message indicating the first n mrna_ids in the genome.

        If no mrnas are present, returns a message to that effect. If fewer than n
        mrnas are loaded, returns the mrna_ids of those mrnas."""
        mrnas_list = []
        while len(mrnas_list) < number:
            for seq in self.seqs:
                mrnas_list.extend(seq.get_mrna_ids())
        # List may now contain more than 'number' ids, or it may contain zero
        if not mrnas_list:
            return "No mrnas currently in memory.\n"
        if len(mrnas_list) > number:
            mrnas_list = mrnas_list[:number]
        result = "First " + str(len(mrnas_list)) + " mrna ids are: "
        result += format_list_with_strings(mrnas_list)
        return result


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
            if len(args) == 1:
                seq_id = args[0]
                for seq in self.seqs:
                    if seq.header == seq_id:
                        return seq.get_subseq()
            elif len(args) == 3:
                seq_id = args[0]
                start = int(args[1])
                stop = int(args[2])
                for seq in self.seqs:
                    if seq.header == seq_id:
                        return seq.get_subseq(start, stop)
            else:
                return "Usage: barfseq <seq_id> <start_index> <end_index>\n"

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
            if self.filter_mgr.dirty or self.seq_fixer.dirty:
                self.stats_mgr.clear_alt()
                sys.stderr.write("Calculating statistics on genome...\n")
                for seq in self.seqs:
                    # Deep copy seq, apply fixes and filters, then update stats
                    cseq = copy.deepcopy(seq)
                    self.seq_fixer.fix(seq)
                    self.filter_mgr.apply_filters(cseq)
                    self.stats_mgr.update_alt(cseq.stats())
                self.filter_mgr.dirty = False
                self.seq_fixer.dirty = False
            return first_line + self.stats_mgr.summary()

## Utility methods

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

    def contains_mrna(self, mrna_id):
        for seq in self.seqs:
            if seq.contains_mrna(mrna_id):
                return True
        return False

    def contains_gene(self, gene_id):
        for seq in self.seqs:
            if seq.contains_gene(gene_id):
                return True
        return False

    def contains_seq(self, seq_id):
        for seq in self.seqs:
            if seq.header == seq_id:
                return True
        return False

    def can_write_to_path(self, path):
        return not os.path.exists(path)


## Utility functions
def format_list_with_strings(entries):
    if len(entries) == 0:
        return ""
    result = entries[0]
    if len(entries) > 2:
        for entry in entries[1:-1]:
            result += ", " + entry
    if len(entries) > 1:
        result += ", " + entries[-1]
    result += "\n"
    return result

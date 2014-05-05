#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import subprocess
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
        self.removed_seqs = []
        self.annot = Annotator()
        self.filter_mgr = FilterManager()
        self.stats_mgr = StatsManager()
        self.seq_fixer = SeqFixer()

    def genome_is_loaded(self):
        for seq in self.seqs:
            if seq.genes:
                return True
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
            gff = open(line+'/genome.gff', 'w')
            removed_gff = open(line+'/genome.removed.gff', 'w')
            tbl = open(line+'/genome.tbl', 'w')
            fasta = open(line+'/genome.fasta', 'w')
            mrna_fasta = open(line+'/genome.mrna.fasta', 'w')
            cds_fasta = open(line+'/genome.cds.fasta', 'w')
            protein_fasta = open(line+'/genome.proteins.fasta', 'w')
            stats_file = open(line+'/genome.stats', 'w')


            # Deep copy each seq, apply fixes and filters, write
            sys.stderr.write("Writing gff, tbl and fasta...\n")
            number_of_gagflags = 0
            first_line = "Number of sequences:   " + str(len(self.seqs)) + "\n"
            update_alt = False
            if self.filter_mgr.dirty or self.seq_fixer.dirty:
                self.stats_mgr.clear_alt()
                update_alt = True
            for seq in self.removed_seqs:
                removed_gff.write(seq.to_gff())
            for seq in self.seqs:
                cseq = copy.deepcopy(seq)
                self.apply_filters_n_fixes(cseq)
                gff.write(cseq.to_gff())
                removed_gff.write(cseq.removed_to_gff())
                tbl.write(cseq.to_tbl())
                mrna_fasta.write(cseq.to_mrna_fasta())
                cds_fasta.write(cseq.to_cds_fasta())
                protein_fasta.write(cseq.to_protein_fasta())
                fasta.write(cseq.to_fasta())
                if update_alt:
                    self.stats_mgr.update_alt(cseq.stats())
                number_of_gagflags += cseq.number_of_gagflags()

            last_line = "(" + str(number_of_gagflags) + " features flagged)\n"
            stats_file.write(first_line + self.stats_mgr.summary() + last_line)
            self.filter_mgr.dirty = False
            self.seq_fixer.dirty = False

            # Close files
            gff.close()
            tbl.close()
            fasta.close()
            mrna_fasta.close()
            cds_fasta.close()
            protein_fasta.close()
            stats_file.close()

            return "Genome written to " + line
        
    def load_folder(self, line):
        if not line:
            line = "."
        fastapath = line + '/genome.fasta'
        gffpath = line + '/genome.gff'

        # Verify files
        if not os.path.isfile(fastapath):
            sys.stderr.write("Failed to find " + fastapath + ". No genome was loaded.")
            return
        if not os.path.isfile(gffpath):
            sys.stderr.write("Failed to find " + gffpath + ". No genome was loaded.")
            return

        # Read the fasta
        sys.stderr.write("Reading fasta...\n")
        self.read_fasta(fastapath)
        sys.stderr.write("Done.\n")

        # Read the gff
        sys.stderr.write("Reading gff...\n")
        self.read_gff(gffpath)
        sys.stderr.write("Done.\n")

        # Clear stats; read in new stats
        self.stats_mgr.clear_all()
        for seq in self.seqs:
            self.stats_mgr.update_ref(seq.stats())

    def add_annotations_from_list(self, anno_list):
        for seq in self.seqs:
            seq.add_annotations_from_list(anno_list)

    def remove_from_file(self, filename):
        if not os.path.isfile(filename):
            sys.stderr.write("Error: " + filename + " is not a file. Nothing removed.\n")
            return
        remove_list = self.read_remove_list(filename)
        self.remove_from_list(remove_list)
        self.filter_mgr.dirty = True

    def trim_from_file(self, filename):
        if not os.path.isfile(filename):
            sys.stderr.write("Error: " + filename + " is not a file. Nothing trimmed.\n")
            return
        trimlist = self.read_bed_file(open(filename, 'rb'))
        if not trimlist:
            sys.stderr.write("Failed to read .bed file; nothing trimmed.\n")
            return
        else:
            self.trim_from_list(trimlist)
            # TODO filter mgr? or seq fixer? or what?
            self.filter_mgr.dirty = True

    def annotate_from_file(self, filename):
        if not os.path.isfile(filename):
            sys.stderr.write("Error: " + filename + " is not a file. Nothing removed.\n")
            return
        annos = self.read_annotation_file(open(filename, 'rb'))
        if not annos:
            sys.stderr.write("Failed to read annotations from " + filename + "; no annotations added.\n")
            return
        else:
            sys.stderr.write("Adding annotations to genome ...\n")
            self.add_annotations_from_list(annos)
            # TODO should make something dirty? we changed the genome...
            self.filter_mgr.dirty = True
            sys.stderr.write("...done\n")

    def trim_from_list(self, trimlist):
        for seq in self.seqs:
            # Trim the ends first
            for entry in trimlist:
                if entry[0] == seq.header and entry[2] == len(seq.bases):
                    seq.trim_region(entry[1], entry[2])
                    sys.stderr.write("Trimmed " + entry[0] + " from ")
                    sys.stderr.write(str(entry[1]) + " to " + str(entry[2]) + "\n")
            # Now trim the beginnings
            for entry in trimlist:
                if entry[0] == seq.header and entry[1] == 1:
                    seq.trim_region(entry[1], entry[2])
                    sys.stderr.write("Trimmed " + entry[0] + " from ")
                    sys.stderr.write(str(entry[1]) + " to " + str(entry[2]) + "\n")

    def set_filter_arg(self, filter_name, val):
        self.filter_mgr.set_filter_arg(filter_name, val)

    def get_filter_arg(self, filter_name):
        return self.filter_mgr.get_filter_arg(filter_name)

    def set_filter_remove(self, filter_name, remove):
        self.filter_mgr.set_filter_remove(filter_name, remove)
        
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

    def read_remove_list(self, line):
        remove_list = []
        with open(line, 'rb') as infile:
            for line in infile:
                remove_list.append(line.strip()) 
        return remove_list

    def read_bed_file(self, io_buffer):
        trimlist = []
        for line in io_buffer:
            splitline = line.strip().split('\t')
            if len(splitline) != 3:
                return []
            else:
                try:
                    entry = [splitline[0], int(splitline[1]), int(splitline[2])]
                except ValueError:
                    sys.stderr.write("Error reading .bed file. Non-integer value ")
                    sys.sdterr.write("in column 2 or 3. Here is the line:\n")
                    sys.stderr.write(line)
                    return []
                trimlist.append(entry)
        return trimlist

    def read_annotation_file(self, io_buffer):
        annos = []
        for line in io_buffer:
            splitline = line.strip().split('\t')
            if len(splitline) != 3:
                return []
            else:
                annos.append(splitline)
        return annos



## Apply filters n fixes

    def apply_filters_n_fixes(self, seq):
        self.seq_fixer.fix(seq)
        self.filter_mgr.apply_filters(seq)
        seq.remove_empty_mrnas()
        seq.remove_empty_genes()
        


## Output info to console

    def barf_gene_gff(self, line):
        if not self.seqs:
            return self.no_genome_message
        else:
            for seq in self.seqs:
                if seq.contains_gene(line):
                    cseq = copy.deepcopy(seq)
                    self.apply_filters_n_fixes(cseq)
                    return cseq.gene_to_gff(line)

    def barf_seq(self, line):
        if not self.seqs:
            return self.no_genome_message
        else:
            args = line.split(' ')
            if len(args) == 1:
                seq_id = args[0]
                for seq in self.seqs:
                    if seq.header == seq_id:
                        cseq = copy.deepcopy(seq)
                        self.apply_filters_n_fixes(cseq)
                        return cseq.get_subseq()
            elif len(args) == 3:
                seq_id = args[0]
                start = int(args[1])
                stop = int(args[2])
                for seq in self.seqs:
                    if seq.header == seq_id:
                        cseq = copy.deepcopy(seq)
                        self.apply_filters_n_fixes(cseq)
                        return cseq.get_subseq(start, stop)
            else:
                return "Usage: barfseq <seq_id> <start_index> <end_index>\n"

    def barf_cds_seq(self, line):
        if not self.seqs:
            return self.no_genome_message
        else:
            name = line
            for seq in self.seqs:
                if seq.contains_mrna(name):
                    cseq = copy.deepcopy(seq)
                    self.apply_filters_n_fixes(cseq)
                    return cseq.extract_cds_seq(name)
            return "Error: Couldn't find mRNA.\n"

    def cds_to_gff(self, line):
        if not self.seqs:
            return self.no_genome_message
        else:
            name = line
            for seq in self.seqs:
                if seq.contains_mrna(name):
                    cseq = copy.deepcopy(seq)
                    self.apply_filters_n_fixes(cseq)
                    return cseq.cds_to_gff(name)
            return "Error: Couldn't find mRNA.\n"

    def cds_to_tbl(self, line):
        if not self.seqs:
            return self.no_genome_message
        else:
            name = line
            for seq in self.seqs:
                if seq.contains_mrna(name):
                    cseq = copy.deepcopy(seq)
                    self.apply_filters_n_fixes(cseq)
                    return cseq.cds_to_tbl(name)
            return "Error: Couldn't find mRNA.\n"

    def barf_gene_tbl(self, line):
        if not self.seqs:
            return self.no_genome_message
        else:
            output = ">Feature SeqId\n"
            for seq in self.seqs:
                if seq.contains_gene(line):
                    cseq = copy.deepcopy(seq)
                    self.apply_filters_n_fixes(cseq)
                    output += cseq.gene_to_tbl(line)
            return output

    def stats(self):
        if not self.seqs:
            return self.no_genome_message
        else:
            number_of_gagflags = 0
            first_line = "Number of sequences:   " + str(len(self.seqs)) + "\n"
            if self.filter_mgr.dirty or self.seq_fixer.dirty:
                self.stats_mgr.clear_alt()
                sys.stderr.write("Calculating statistics on genome...\n")
                for seq in self.seqs:
                    # Deep copy seq, apply fixes and filters, then update stats
                    cseq = copy.deepcopy(seq)
                    self.apply_filters_n_fixes(cseq)
                    self.stats_mgr.update_alt(cseq.stats())
                    number_of_gagflags += cseq.number_of_gagflags()
                self.filter_mgr.dirty = False
                self.seq_fixer.dirty = False
            last_line = "(" + str(number_of_gagflags) + " features flagged)\n"
            return first_line + self.stats_mgr.summary() + last_line

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

    def remove_from_list(self, bad_list):
        # First remove any seqs on the list
        to_remove = []
        for seq in self.seqs:
            if seq.header in bad_list:
                to_remove.append(seq)
        if to_remove:
            for seq in to_remove:
                self.seqs.remove(seq)
                # TODO this is where the whole deepcopy scheme starts to unravel
                # we're handling removed seqs differently from other removed features.
                sys.stderr.write("Warning: removing seq " + seq.header + ".\n")
                sys.stderr.write("You must reload genome to get this sequence back.\n")
            self.removed_seqs.extend(to_remove)
        # Now pass the list down to each seq
        for seq in self.seqs:
            seq.remove_from_list(bad_list)

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
        if len(path.split()) > 1:
            return False
        else:
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

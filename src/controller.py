#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
from src.fasta_reader import FastaReader
from src.gff_reader import GFFReader
from src.annotator import Annotator
from src.filter_manager import FilterManager
from src.stats_manager import StatsManager

class Controller:

    def __init__(self):
        self.seqs = []
        self.removed_features = []
        self.annot = Annotator()
        self.filter_mgr = FilterManager()
        self.stats_mgr = StatsManager()

    def execute(self, args_dict):
        """At a minimum, write a fasta, gff and tbl to output directory. Optionally do more."""
        # Verify and read fasta file
        fastapath = args_dict["fasta"]
        if not os.path.isfile(fastapath):
            sys.stderr.write("Failed to find " + fastapath + ". No genome was loaded.\n")
            sys.exit()
        sys.stderr.write("Reading fasta...\n")
        self.read_fasta(fastapath)
        sys.stderr.write("Done.\n")

        # Verify and read gff file
        gffpath = args_dict["gff"]
        if not os.path.isfile(gffpath):
            sys.stderr.write("Failed to find " + gffpath + ". No genome was loaded.")
            return
        sys.stderr.write("Reading gff...\n")
        self.read_gff(gffpath)
        sys.stderr.write("Done.\n")

        # Optional annotation step
        if "anno" in args_dict:
            anno_filename = args_dict["anno"]
            self.annotate_from_file(anno_filename)

        # Write fasta, gff and tbl file to gag_output/ folder
        # TODO look for 'out' arg
        out_dir = "gag_output"
        # Create directory, open files
        os.system('mkdir ' + out_dir)
        # Open files
        fasta = open(out_dir+'/genome.fasta', 'w')
        gff = open(out_dir+'/genome.gff', 'w')
        tbl = open(out_dir+'/genome.tbl', 'w')
        sys.stderr.write("Writing gff, tbl and fasta...\n")
        # TODO track # of gagflags?
        # TODO stats file
        for seq in self.seqs:
            fasta.write(seq.to_fasta())
            gff.write(seq.to_gff())
            tbl.write(seq.to_tbl())
        # Close files
        gff.close()
        tbl.close()
        fasta.close()

    def barf_folder(self, line):
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

        # Now write stuff
        sys.stderr.write("Writing gff, tbl and fasta...\n")
        number_of_gagflags = 0
        first_line = "Number of sequences:   " + str(len(self.seqs)) + "\n"
        update_alt = False
        self.stats_mgr.clear_alt()
        removed_gff.write("##gff-version 3\n")
        gff.write("##gff-version 3\n")
        for feature in self.removed_features:
            removed_gff.write(feature.to_gff())
        for seq in self.seqs:
            gff.write(seq.to_gff())
            removed_gff.write(seq.removed_to_gff())
            tbl.write(seq.to_tbl())
            mrna_fasta.write(seq.to_mrna_fasta())
            cds_fasta.write(seq.to_cds_fasta())
            protein_fasta.write(seq.to_protein_fasta())
            fasta.write(seq.to_fasta())
            self.stats_mgr.update_alt(seq.stats())
            number_of_gagflags += seq.number_of_gagflags()

        last_line = "(" + str(number_of_gagflags) + " features flagged)\n"
        stats_file.write(first_line + self.stats_mgr.summary() + last_line)

        # Close files
        gff.close()
        tbl.close()
        fasta.close()
        mrna_fasta.close()
        cds_fasta.close()
        protein_fasta.close()
        stats_file.close()
        sys.stderr.write( "Genome written to " + line + "\n")
        
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

        # Remove empty features
        for seq in self.seqs:
            self.remove_empty_features(seq)

        # Clear stats; read in new stats
        sys.stderr.write("Calculating stats...\n")
        self.stats_mgr.clear_all()
        for seq in self.seqs:
            self.stats_mgr.update_ref(seq.stats())
        sys.stderr.write("Done.\n")
    
    def add_annotations_from_list(self, anno_list):
        for seq in self.seqs:
            seq.add_annotations_from_list(anno_list)

    def remove_from_file(self, filename):
        if not os.path.isfile(filename):
            sys.stderr.write("Error: " + filename + " is not a file. Nothing removed.\n")
            return
        remove_list = self.read_remove_list(filename)
        self.remove_from_list(remove_list)

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

    def annotate_from_file(self, filename):
        if not os.path.isfile(filename):
            sys.stderr.write("Error: " + filename + " is not a file. Nothing annotated.\n")
            return
        annos = self.read_annotation_file(open(filename, 'rb'))
        if not annos:
            sys.stderr.write("Failed to read annotations from " + filename + "; no annotations added.\n")
            return
        else:
            sys.stderr.write("Adding annotations to genome ...\n")
            self.add_annotations_from_list(annos)
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
            self.remove_empty_features(seq)

    def get_filter_arg(self, filter_name):
        return self.filter_mgr.get_filter_arg(filter_name)
        
    def apply_filter(self, filter_name, val, filter_mode):
        for seq in self.seqs:
            self.filter_mgr.apply_filter(filter_name, val, filter_mode, seq)
            self.remove_empty_features(seq)

    def fix_terminal_ns(self):
        for seq in self.seqs:
            seq.remove_terminal_ns()
            self.remove_empty_features(seq)
        return "Terminal Ns fixed."

    def fix_start_stop_codons(self):
        for seq in self.seqs:
            seq.create_starts_and_stops()
        return "Verified and created start/stop codons."

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


## Clean up

    def remove_empty_features(self, seq):
        """Removes any empty mRNAs or genes from a seq and adds them to self.removed_features."""
        self.removed_features.extend(seq.remove_empty_mrnas())
        self.removed_features.extend(seq.remove_empty_genes())
        
    def stats(self):
        if not self.seqs:
            return self.no_genome_message
        else:
            number_of_gagflags = 0
            # TODO have stats mgr handle "number of sequences"
            first_line = "Number of sequences:   " + str(len(self.seqs)) + "\n"
            sys.stderr.write("Calculating statistics on genome...\n")
            self.stats_mgr.clear_alt()
            for seq in self.seqs:
                self.stats_mgr.update_alt(seq.stats())
                number_of_gagflags += seq.number_of_gagflags()
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
                sys.stderr.write("Warning: removing seq " + seq.header + ".\n")
                sys.stderr.write("You must reload genome to get this sequence back.\n")
            self.removed_features.extend(to_remove)
        # Now pass the list down to each seq
        for seq in self.seqs:
            removed_from_seq = seq.remove_from_list(bad_list)
            self.removed_features.extend(removed_from_seq)

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

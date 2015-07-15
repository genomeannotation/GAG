#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import sys
from src.controller import Controller

def usage():
    sys.stderr.write("usage:\tpython gag.py <fasta=fasta_file> <gff=gff_file> ")
    sys.stderr.write("\n\t[anno=<annotation file>] [trim=<bed file>] ")
    sys.stderr.write("\n\t[out=<output directory>] [fix_start_stop=true/false] ")
    sys.stderr.write("\n\t[fix_terminal_ns=true/false] ")
    sys.stderr.write("\n\t[remove_cds_shorter_than=<length>] [remove_cds_longer_than=<length>] ")
    sys.stderr.write("\n\t[remove_exons_shorter_than=<length>] [remove_exons_longer_than=<length>] ")
    sys.stderr.write("\n\t[remove_introns_shorter_than=<length>] [remove_introns_longer_than=<length>] ")
    sys.stderr.write("\n\t[remove_genes_shorter_than=<length>] [remove_genes_longer_than=<length>] ")
    sys.stderr.write("\n\t[flag_cds_shorter_than=<length>] [flag_cds_longer_than=<length>] ")
    sys.stderr.write("\n\t[flag_exons_shorter_than=<length>] [flag_exons_longer_than=<length>] ")
    sys.stderr.write("\n\t[flag_introns_shorter_than=<length>] [flag_introns_longer_than=<length>] ")
    sys.stderr.write("\n\t[flag_genes_shorter_than=<length>] [flag_genes_longer_than=<length>] ")
    sys.stderr.write("\n")
    sys.exit()

def parse_args():
    if len(sys.argv) == 1 or sys.argv[1] == "help":
        usage()
    have_fasta = False
    have_gff = False
    for arg in sys.argv:
        if "fasta" in arg:
            have_fasta = True
            continue
        elif "gff" in arg:
            have_gff = True
    if not have_fasta or not have_gff:
        sys.stderr.write("Error: must provide a fasta and gff file.\n")
        usage()
    args = {}
    for arg in sys.argv:
        if "=" in arg:
            split_arg = arg.split("=")
            args[split_arg[0]] = split_arg[1]
    return args

def main():
    args = parse_args()
    controller = Controller()
    controller.execute(args)

if __name__ == '__main__':
    main()

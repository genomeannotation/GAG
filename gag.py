#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import sys
from src.controller import Controller

def usage():
    sys.stderr.write("usage: python gag.py <fasta=fasta_file> <gff=gff_file> ")
    sys.stderr.write("[anno=<annotation file>] [trim=<bed file>] ")
    sys.stderr.write("[out=<output directory>] [fix_start_stop=true/false]")
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

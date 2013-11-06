#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import csv
from src.feature_classes import *
from src.gff import GFF
from src.bed import Bed

gff_filename = "demo/demo.gff"
bed_filename = "demo/demo.bed"

gff = GFF()
bed = Bed()

with open(gff_filename, 'rb') as gfffile:
    gffreader = csv.reader(gfffile, delimiter='\t')
    gff.read_file(gffreader)

with open(bed_filename, 'rb') as bedfile:
    bedreader = csv.reader(bedfile, delimiter='\t')
    bed.read_file(bedreader)


# TODO parse bed and apply it...

print("Trimming first gene (on sctg_1)...")
print("Here is the gff output:\n")

for gene in gff.genes:
    if gene.seq_name == 'sctg_1':
        gene.trim([201, 299])
    print(gene.to_gff())

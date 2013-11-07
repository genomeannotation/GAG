#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import csv
from src.gene import Gene
from src.mrna import MRNA
from src.gene_part import GenePart, CDS, Exon
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


print("here is the gff I read...")
for gene in gff.genes:
    print(gene.to_gff())

print("\n")
print("here is the bed...")
print(str(bed))

print("okay, last gene's first mrna:")
print(str(gff.genes[len(gff.genes)-1].mrnas[0]))

gff.apply_bed(bed)
gff.remove_empty_genes()

print("okay, last gene's first mrna:")
print(str(gff.genes[len(gff.genes)-1].mrnas[0]))
print("\n")
print("and here is the gff after applying changes")
for gene in gff.genes:
    print(gene.to_gff())

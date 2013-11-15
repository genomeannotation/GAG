#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import csv
from src.fasta import Fasta
from src.gene import Gene
from src.mrna import MRNA
from src.gene_part import GenePart, CDS, Exon
from src.gff import GFF
from src.bed import Bed
from src.feature_tbl_entry import FeatureTblEntry
from src.annotator import Annotator

#gff_filename = "demo/demo.gff"
bed_filename = "demo/demo.bed"
#fasta_filename = "demo/demo.fasta"
#trinotate_filename = "demo/demo.trinotate"

gff_filename = "real_files/real_output.gff"
fasta_filename = "real_files/454ScaffoldContigs.fna"
trinotate_filename = "real_files/maker.xls"

gff = GFF()
bed = Bed()
fasta = Fasta()

fasta.read_file(fasta_filename)

with open(gff_filename, 'rb') as gfffile:
    gffreader = csv.reader(gfffile, delimiter='\t')
    gff.read_file(gffreader)

with open(bed_filename, 'rb') as bedfile:
    bedreader = csv.reader(bedfile, delimiter='\t')
    bed.read_file(bedreader)

print("********BEFORE\n*******************")
#print("Fasta:\n")
#print(fasta.write_string())
#print("\nGFF:\n")
#for gene in gff.genes:
#    sys.stdout.write(gene.to_gff()) 

print("\n\n********TABLE WRITING\n*******************")
outFile = open('out.tbl', 'w')
outFile.write('>Feature SeqId\n')

annot = Annotator()
annot.read_from_file(trinotate_filename)
for seq in fasta.entries:
    outFile.write('>Feature '+seq[0]+'\n')
    outFile.write('1\t'+str(len(seq[1]))+'\tREFERENCE\n\t\t\tPBARC\t12345\n')
    for gene in gff.genes:
        if gene.seq_name != seq[0]:
            continue

        entries = gene.to_tbl_entries()
        for entry in entries:    
            if entry.type == 'gene':
                annot.annotate_gene(entry)
            elif entry.type == 'CDS':
                annot.annotate_cds(entry)
            elif entry.type == 'mRNA':
                annot.annotate_mrna(entry)
            outFile.write(entry.write_to_string()+'\n')
            #print(entry.write_to_string())
outFile.close()

exit(0)

fasta.apply_bed(bed)
gff.apply_bed(bed)
gff.remove_empty_genes()

print("\n********AFTER\n*******************")
print("Fasta:\n")
print(fasta.write_string())
print("\nGFF:\n")
for gene in gff.genes:
    sys.stdout.write(gene.to_gff())

print("\n\n********TABLE WRITING\n*******************")
annot = Annotator()
annot.read_from_file(trinotate_filename)
for gene in gff.genes:
    entries = gene.to_tbl_entries()
    for entry in entries:    
        if entry.type == 'gene':
            annot.annotate_gene(entry)
        elif entry.type == 'CDS':
            annot.annotate_cds(entry)
        elif entry.type == 'mRNA':
            annot.annotate_mrna(entry)
        print(entry.write_to_string())

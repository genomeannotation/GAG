#!/usr/bin/env python

class CDSLengthRangeFilter:

    def __init__(self, min_length = 0, max_length=0):
        self.min_length = min_length
        self.max_length = max_length
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if not mrna.cds:
                    gene.add_annotation("invalidated_mrna", mrna.identifier+" didn't pass cds_length_range - CDS doesn't exist")
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on
                elif mrna.cds.length() < self.min_length:
                    mrna.cds.add_annotation("invalidated", "didn't pass cds_length_range with CDS shorter than "+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on
                elif self.max_length > 0 and mrna.cds.length() > self.max_length:
                    mrna.cds.add_annotation("invalidated", "didn't pass cds_length_range with CDS longer than "+str(self.max_length))
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on
                    
class ExonLengthRangeFilter:

    def __init__(self, min_length = 0, max_length=0):
        self.min_length = min_length
        self.max_length = max_length
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if not mrna.exon:
                    gene.add_annotation("invalidated_mrna", mrna.identifier+" didn't pass exon_length_range - Exon doesn't exist")
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on
                elif mrna.get_shortest_exon() < self.min_length:
                    mrna.exon.add_annotation("invalidated", "didn't pass exon_length_range with Exon shorter than "+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on
                elif self.max_length > 0 and mrna.get_longest_exon() > self.max_length:
                    mrna.exon.add_annotation("invalidated", "didn't pass exon_length_range with Exon longer than "+str(self.max_length))
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on

class IntronLengthRangeFilter:

    def __init__(self, min_length = 0, max_length=0):
        self.min_length = min_length
        self.max_length = max_length
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if not mrna.exon:
                    gene.add_annotation("invalidated_mrna", mrna.identifier+" didn't pass intron_length_range - Exon doesn't exist")
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on
                elif mrna.get_shortest_intron() < self.min_length:
                    mrna.exon.add_annotation("invalidated", "didn't pass intron_length_range with intron shorter than "+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the intron lives on
                elif self.max_length > 0 and mrna.get_longest_intron() > self.max_length:
                    mrna.exon.add_annotation("invalidated", "didn't pass intron_length_range with intron longer than "+str(self.max_length))
                    mrna.death_flagged = True # Destroy the mRNA that the intron lives on

#!/usr/bin/env python

class CDSLengthRangeFilter:

    def __init__(self, min_length = 0, max_length=0):
        self.min_length = min_length
        self.max_length = max_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if not mrna.cds:
                    gene.add_annotation('cds_length_range_invalidated_mrna', mrna.identifier+" didn't pass - CDS doesn't exist")
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on?
                elif mrna.cds.length() < self.min_length:
                    mrna.cds.add_annotation('cds_length_range_invalidated', "didn't pass with CDS shorter than "+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on?
                elif self.max_length > 0 and mrna.cds.length() > self.max_length:
                    mrna.cds.add_annotation('cds_length_range_invalidated', "didn't pass with CDS longer than "+str(self.max_length))
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on?
            if self.remove:
                gene.mrnas = [mrna for mrna in gene.mrnas if not mrna.death_flagged]
            for mrna in gene.mrnas:
                mrna.death_flagged = False
  
                    
class ExonLengthRangeFilter:

    def __init__(self, min_length = 0, max_length=0):
        self.min_length = min_length
        self.max_length = max_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if not mrna.exon:
                    gene.add_annotation("exon_length_range_invalidated", mrna.identifier+" didn't pass - Exon doesn't exist")
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on?
                elif mrna.get_shortest_exon() < self.min_length:
                    mrna.exon.add_annotation("exon_length_range_invalidated", "didn't pass with Exon shorter than "+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on?
                elif self.max_length > 0 and mrna.get_longest_exon() > self.max_length:
                    mrna.exon.add_annotation("exon_length_range_invalidated", "didn't pass with Exon longer than "+str(self.max_length))
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on?
            if self.remove:
                gene.mrnas = [mrna for mrna in gene.mrnas if not mrna.death_flagged]
            for mrna in gene.mrnas:
                mrna.death_flagged = False

class IntronLengthRangeFilter:

    def __init__(self, min_length = 0, max_length=0):
        self.min_length = min_length
        self.max_length = max_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if not mrna.exon:
                    gene.add_annotation("intron_length_range_invalidated", mrna.identifier+" didn't pass - Exon doesn't exist")
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on?
                elif mrna.get_shortest_intron() < self.min_length:
                    mrna.exon.add_annotation("intron_length_range_invalidated", "didn't pass with intron shorter than "+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the intron lives on?
                elif self.max_length > 0 and mrna.get_longest_intron() > self.max_length:
                    mrna.exon.add_annotation("intron_length_range_invalidated", "didn't pass with intron longer than "+str(self.max_length))
                    mrna.death_flagged = True # Destroy the mRNA that the intron lives on?
            if self.remove:
                gene.mrnas = [mrna for mrna in gene.mrnas if not mrna.death_flagged]
            for mrna in gene.mrnas:
                mrna.death_flagged = False
                    
class GeneLengthRangeFilter:

    def __init__(self, min_length = 0, max_length=0):
        self.min_length = min_length
        self.max_length = max_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            if gene.length() < self.min_length:
                gene.add_annotation("gene_length_range_invalidated", "didn't pass with gene shorter than "+str(self.min_length))
                gene.death_flagged = True # Destroy the gene?
            elif self.max_length > 0 and gene.length() > self.max_length:
                gene.add_annotation("gene_length_range_invalidated", "didn't pass with gene longer than "+str(self.max_length))
                gene.death_flagged = True # Destroy the gene?
        if self.remove:
            seq.genes = [gene for gene in seq.genes if not gene.death_flagged]
        for gene in seq.genes:
            gene.death_flagged = False

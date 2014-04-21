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
                    gene.add_annotation('gag_flag', "cds_min_length:"+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on?
                elif mrna.cds.length() < self.min_length:
                    mrna.cds.add_annotation('gag_flag', "cds_min_length:"+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on?
                elif self.max_length > 0 and mrna.cds.length() > self.max_length:
                    mrna.cds.add_annotation('gag_flag', "cds_max_length:"+str(self.max_length))
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on?
            if self.remove:
                gene.mrnas = [mrna for mrna in gene.mrnas if not mrna.death_flagged]
            for mrna in gene.mrnas:
                mrna.death_flagged = False
            # Mark empty genes for removal
            if not gene.mrnas:
                gene.death_flagged = True
        # Remove empty genes
        seq.genes = [g for g in seq.genes if not g.death_flagged]
  
                    
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
                    gene.add_annotation("gag_flag", "exon_min_length:"+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on?
                elif mrna.get_shortest_exon() < self.min_length:
                    mrna.exon.add_annotation("gag_flag", "exon_min_length:"+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on?
                elif self.max_length > 0 and mrna.get_longest_exon() > self.max_length:
                    mrna.exon.add_annotation("gag_flag", "exon_max_length:"+str(self.max_length))
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on?
            if self.remove:
                gene.mrnas = [mrna for mrna in gene.mrnas if not mrna.death_flagged]
            for mrna in gene.mrnas:
                mrna.death_flagged = False
            # Mark empty genes for removal
            if not gene.mrnas:
                gene.death_flagged = True
        # Remove empty genes
        seq.genes = [g for g in seq.genes if not g.death_flagged]

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
                    gene.add_annotation("gag_flag", "intron_min_length:"+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on?
                elif mrna.get_shortest_intron() < self.min_length:
                    mrna.exon.add_annotation("gag_flag", "intron_min_length:"+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the intron lives on?
                elif self.max_length > 0 and mrna.get_longest_intron() > self.max_length:
                    mrna.exon.add_annotation("gag_flag", "intron_max_length:"+str(self.max_length))
                    mrna.death_flagged = True # Destroy the mRNA that the intron lives on?
            if self.remove:
                gene.mrnas = [mrna for mrna in gene.mrnas if not mrna.death_flagged]
            for mrna in gene.mrnas:
                mrna.death_flagged = False
            # Mark empty genes for removal
            if not gene.mrnas:
                gene.death_flagged = True
        # Remove empty genes
        seq.genes = [g for g in seq.genes if not g.death_flagged]
                    
class GeneLengthRangeFilter:

    def __init__(self, min_length = 0, max_length=0):
        self.min_length = min_length
        self.max_length = max_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            if gene.length() < self.min_length:
                gene.add_annotation("gag_flag", "gene_min_length:"+str(self.min_length))
                gene.death_flagged = True # Destroy the gene?
            elif self.max_length > 0 and gene.length() > self.max_length:
                gene.add_annotation("gag_flag", "gene_max_length:"+str(self.max_length))
                gene.death_flagged = True # Destroy the gene?
        if self.remove:
            seq.genes = [gene for gene in seq.genes if not gene.death_flagged]
        for gene in seq.genes:
            gene.death_flagged = False

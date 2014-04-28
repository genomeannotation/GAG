#!/usr/bin/env python

###################################################################################################

class MinCDSLengthFilter:

    def __init__(self, min_length = 0):
        self.arg = min_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.cds and mrna.cds.length() < self.arg:
                    mrna.cds.add_annotation('gag_flag', "cds_min_length:"+str(self.arg))
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on?
            if self.remove:
                to_remove = [mrna for mrna in gene.mrnas if mrna.death_flagged]
                for mrna in to_remove:
                    gene.remove_mrna(mrna.identifier)
            for mrna in gene.mrnas:
                mrna.death_flagged = False
            # Mark empty genes for removal
            if not gene.mrnas:
                gene.death_flagged = True
        # Remove empty genes
        to_remove = [g for g in seq.genes if g.death_flagged]
        for gene in to_remove:
            seq.remove_gene(gene.identifier)
        
class MaxCDSLengthFilter:

    def __init__(self, max_length=0):
        self.arg = max_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.cds and self.arg > 0 and mrna.cds.length() > self.arg:
                    mrna.cds.add_annotation('gag_flag', "cds_max_length:"+str(self.arg))
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on?
            if self.remove:
                to_remove = [mrna for mrna in gene.mrnas if mrna.death_flagged]
                for mrna in to_remove:
                    gene.remove_mrna(mrna.identifier)
            for mrna in gene.mrnas:
                mrna.death_flagged = False
            # Mark empty genes for removal
            if not gene.mrnas:
                gene.death_flagged = True
        # Remove empty genes
        to_remove = [g for g in seq.genes if g.death_flagged]
        for gene in to_remove:
            seq.remove_gene(gene.identifier)
        
###################################################################################################
                    
class MinExonLengthFilter:

    def __init__(self, min_length = 0):
        self.arg = min_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.exon and mrna.get_shortest_exon() < self.arg:
                    mrna.exon.add_annotation("gag_flag", "exon_min_length:"+str(self.arg))
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on?
            if self.remove:
                to_remove = [mrna for mrna in gene.mrnas if mrna.death_flagged]
                for mrna in to_remove:
                    gene.remove_mrna(mrna.identifier)
            for mrna in gene.mrnas:
                mrna.death_flagged = False
            # Mark empty genes for removal
            if not gene.mrnas:
                gene.death_flagged = True
        # Remove empty genes
        to_remove = [g for g in seq.genes if g.death_flagged]
        for gene in to_remove:
            seq.remove_gene(gene.identifier)
        
class MaxExonLengthFilter:

    def __init__(self, max_length=0):
        self.arg = max_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.exon and self.arg > 0 and mrna.get_longest_exon() > self.arg:
                    mrna.exon.add_annotation("gag_flag", "exon_max_length:"+str(self.arg))
                    mrna.death_flagged = True # Destroy the mRNA that the exon lives on?
            if self.remove:
                to_remove = [mrna for mrna in gene.mrnas if mrna.death_flagged]
                for mrna in to_remove:
                    gene.remove_mrna(mrna.identifier)
            for mrna in gene.mrnas:
                mrna.death_flagged = False
            # Mark empty genes for removal
            if not gene.mrnas:
                gene.death_flagged = True
        # Remove empty genes
        to_remove = [g for g in seq.genes if g.death_flagged]
        for gene in to_remove:
            seq.remove_gene(gene.identifier)

###################################################################################################

class MinIntronLengthFilter:

    def __init__(self, min_length = 0, max_length=0):
        self.arg = min_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.exon and mrna.get_shortest_intron() < self.arg and mrna.get_shortest_intron() != 0:
                    print("\n\n\nmrna.exon and getshortestintron < self.arg...")
                    mrna.exon.add_annotation("gag_flag", "intron_min_length:"+str(self.arg))
                    mrna.death_flagged = True # Destroy the mRNA that the intron lives on?
            if self.remove:
                to_remove = [mrna for mrna in gene.mrnas if mrna.death_flagged]
                print(to_remove)
                for mrna in to_remove:
                    gene.remove_mrna(mrna.identifier)
            for mrna in gene.mrnas:
                mrna.death_flagged = False
            # Mark empty genes for removal
            if not gene.mrnas:
                print("\n\n\nnot gene.mrnas")
                gene.death_flagged = True
        # Remove empty genes
        to_remove = [g for g in seq.genes if g.death_flagged]
        for gene in to_remove:
            seq.remove_gene(gene.identifier)

class MaxIntronLengthFilter:

    def __init__(self, max_length=0):
        self.arg = max_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.exon and self.arg > 0 and mrna.get_longest_intron() > self.arg:
                    mrna.exon.add_annotation("gag_flag", "intron_max_length:"+str(self.arg))
                    mrna.death_flagged = True # Destroy the mRNA that the intron lives on?
            if self.remove:
                to_remove = [mrna for mrna in gene.mrnas if mrna.death_flagged]
                for mrna in to_remove:
                    gene.remove_mrna(mrna.identifier)
            for mrna in gene.mrnas:
                mrna.death_flagged = False
            # Mark empty genes for removal
            if not gene.mrnas:
                gene.death_flagged = True
        # Remove empty genes
        to_remove = [g for g in seq.genes if g.death_flagged]
        for gene in to_remove:
            seq.remove_gene(gene.identifier)

###################################################################################################

class MinGeneLengthFilter:

    def __init__(self, min_length = 0):
        self.arg = min_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            if gene.length() < self.arg:
                gene.add_annotation("gag_flag", "gene_min_length:"+str(self.arg))
                gene.death_flagged = True # Destroy the gene?
        if self.remove:
            seq.genes = [gene for gene in seq.genes if not gene.death_flagged]
        for gene in seq.genes:
            gene.death_flagged = False

class MaxGeneLengthFilter:

    def __init__(self, max_length=0):
        self.arg = max_length
        self.remove = True
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            if self.arg > 0 and gene.length() > self.arg:
                gene.add_annotation("gag_flag", "gene_max_length:"+str(self.arg))
                gene.death_flagged = True # Destroy the gene?
        if self.remove:
            seq.genes = [gene for gene in seq.genes if not gene.death_flagged]
        for gene in seq.genes:
            gene.death_flagged = False

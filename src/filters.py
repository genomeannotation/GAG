#!/usr/bin/env python

class CDSLengthRangeFilter:

    def __init__(self, min_length = 0, max_length=0):
        self.min_length = min_length
        self.max_length = max_length
        return
        
    def apply(self, seq):
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.cds != None and mrna.cds.length() < self.min_length:
                    mrna.cds.add_annotation("invalidated", "didn't pass cds_length_range with CDS shorter than "+str(self.min_length))
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on
                elif mrna.cds != None and self.max_length > 0 and mrna.cds.length() > self.max_length:
                    mrna.cds.add_annotation("invalidated", "didn't pass cds_length_range with CDS longer than "+str(self.max_length))
                    mrna.death_flagged = True # Destroy the mRNA that the cds lives on

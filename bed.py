#!/usr/bin/env python

# Bed contains a dict called 'entries' whose keys are chromosomes/contigs
# and whose values are lists of coordinates
# Example: awesome_entry = {'sctg_0002_0347': [1, 4307]}
class Bed:


    def __init__(self, entries=None):
        if entries == None:
            self.entries = {}
        else:
            self.entries = entries

    def add_entry(self, key, val):
        self.entries[key] = val

    def contains(self, chrom):
        return chrom in self.entries

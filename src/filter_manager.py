#!/usr/bin/env python

import ast

# All the filters
from src.filters import *

class FilterManager:

    def __init__(self):
        # Build filters
        self.filters = dict()
        self.filters['cds_shorter_than'] = MinCDSLengthFilter()
        self.filters['cds_longer_than'] = MaxCDSLengthFilter()
        self.filters['exon_shorter_than'] = MinExonLengthFilter()
        self.filters['exon_longer_than'] = MaxExonLengthFilter()
        self.filters['intron_shorter_than'] = MinIntronLengthFilter()
        self.filters['intron_longer_than'] = MaxIntronLengthFilter()
        self.filters['gene_shorter_than'] = MinGeneLengthFilter()
        self.filters['gene_longer_than'] = MaxGeneLengthFilter()
        
        # Starts out dirty
        self.dirty = False

    def apply_filters(self, seq):
        for filt in self.filters.values():
            filt.apply(seq)
    
    def set_filter_arg(self, filter_name, val):
        val = ast.literal_eval(val)
        if self.filters[filter_name].arg != val:
            self.dirty = True
        self.filters[filter_name].arg = val
    
    def get_filter_arg(self, filter_name):
        return self.filters[filter_name].arg
    
    def set_filter_remove(self, filter_name, remove):
        if self.filters[filter_name].remove != remove:
            self.dirty = True
        self.filters[filter_name].remove = remove
   

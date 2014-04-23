#!/usr/bin/env python

import ast

# All the filters
from src.filters import *

class FilterManager:

    def __init__(self):
        # Build filters
        self.filters = dict()
        self.filters['min_cds_length'] = MinCDSLengthFilter()
        self.filters['max_cds_length'] = MaxCDSLengthFilter()
        self.filters['min_exon_length'] = MinExonLengthFilter()
        self.filters['max_exon_length'] = MaxExonLengthFilter()
        self.filters['min_intron_length'] = MinIntronLengthFilter()
        self.filters['max_intron_length'] = MaxIntronLengthFilter()
        self.filters['min_gene_length'] = MinGeneLengthFilter()
        self.filters['max_gene_length'] = MaxGeneLengthFilter()
        
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
        self.filters[filter_name].remove = remove
   

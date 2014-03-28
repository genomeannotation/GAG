#!/usr/bin/env python

import ast

# All the filters
from src.min_cds_length_filter import MinCDSLengthFilter

class FilterManager:

    def __init__(self):
        self.filters = dict()
        self.filters['min_cds_length'] = MinCDSLengthFilter()
        
    def apply_filters(self, seq):
        for filt in self.filters:
            filt.apply(seq)
    
    def set_filter_arg(self, filter_name, filter_arg, val):
        setattr(self.filters[filter_name], filter_arg, ast.literal_eval(val))
    
    def get_filter_arg(self, filter_name, filter_arg):
        return getattr(self.filters[filter_name], filter_arg)
   

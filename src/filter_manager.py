#!/usr/bin/env python

# All the filters
from src.min_cds_length_filter import MinCDSLengthFilter

class FilterManager:

    def __init__(self):
        self.filters = []
        self.filters.append(MinCDSLengthFilter())
        
    def apply_filters(self, seq):
        for filt in filters:
            filt.apply(seq)

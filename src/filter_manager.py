#!/usr/bin/env python

class FilterManager:

    def __init__(self):
        self.filters = []
        return
        
    def add_filter(self, filt):
        self.filters.append(filt)
        
    def apply_filters(self, seq):
        for filt in filters:
            filt.apply(seq)

#!/usr/bin/env python

import ast

# All the filters
from src.filters import CDSLengthRangeFilter

class FilterManager:

    def __init__(self):
        # Build filters
        self.filters = dict()
        self.filters['cds_length_range'] = CDSLengthRangeFilter()
        
        # Build args
        self.filter_args = dict()
        for filt_name, filt in self.filters.items():
            self.filter_args[filt_name] = [attr for attr in dir(filt) if not callable(getattr(filt, attr)) and not attr.startswith("__")]
        
        # Starts out dirty
        self.dirty = True

    def apply_filters(self, seq):
        for filt in self.filters.values():
            filt.apply(seq)
    
    def set_filter_arg(self, filter_name, filter_arg, val):
        val = ast.literal_eval(val)
        if self.get_filter_arg(filter_name, filter_arg) != val:
            self.dirty = True
        setattr(self.filters[filter_name], filter_arg, val)
    
    def get_filter_arg(self, filter_name, filter_arg):
        return getattr(self.filters[filter_name], filter_arg)
        
    def get_filter_help(self, filter_name = ""):
        # If there is no filter name, supply all of the filters' args
        if filter_name == "":
            args_help = ''
            for filt_name, args in self.filter_args.items():
                args_help += filt_name+": " + ", ".join(args)
            return args_help
        elif filter_name in self.filters:
            return filter_name+": " + ", ".join(self.filter_args[filter_name])
        else:
            return "No such filter: "+filter_name
   

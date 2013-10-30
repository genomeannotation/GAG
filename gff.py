#!/usr/bin/env python

def validate_entry(key, val):
    if isinstance(val, list) and isinstance(val[0], int,) \
    and isinstance(val[1], int) and not isinstance(key, list):
        return True
    else:
        return False

class GFF:


    def __init__(self, entries=None):
        if entries == None:
            self.entries = {}
        else:
            self.entries = entries


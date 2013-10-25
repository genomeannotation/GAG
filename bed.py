#!/usr/bin/env python

# Bed contains a dict called 'entries' whose keys are chromosomes/contigs
# and whose values are lists of coordinates
# Example: awesome_entry = {'sctg_0002_0347': [1, 4307]}

def validate_entry(key, val):
    if isinstance(val, list) and isinstance(val[0], int,) \
    and isinstance(val[1], int) and not isinstance(key, list):
        return True
    else:
        return False

class Bed:


    def __init__(self, entries=None):
        if entries == None:
            self.entries = {}
        else:
            self.entries = entries

    def add_entry(self, key, val):
        if validate_entry(key, val):
            self.entries[key] = val
        else:
            raise TypeError(key)

    def contains(self, seq_id):
        return seq_id in self.entries

    def get_coordinates(self, seq_id):
        if seq_id in self.entries:
            return self.entries[seq_id]
        else:
            return None

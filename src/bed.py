#!/usr/bin/env python

# Bed contains a dict called 'entries' whose keys are chromosomes/contigs
# and whose values are lists of coordinates
# Example: awesome_entry = {'sctg_0002_0347': [1, 4307]}

import sys

def validate_entry(key, val):
    if isinstance(val, list) and isinstance(val[0], int,) \
    and isinstance(val[1], int) and not isinstance(key, list):
        return True
    else:
        return False

class Bed:


    def __init__(self, entries=None):
        self.entries = {}
        if entries is not None:
            self.entries.update(entries)

    def __str__(self):
        result = "Bed containing "
        result += str(len(self.entries))
        result += " entries"
        return result

    def add_entry(self, key, val):
        if validate_entry(key, val):
            self.entries[key] = val
        else:
            raise TypeError(key, val)

    def contains(self, seq_id):
        return seq_id in self.entries

    def get_coordinates(self, seq_id):
        if seq_id in self.entries:
            return self.entries[seq_id]
        else:
            return None

    def process_line(self, line):
        if len(line) < 3:
            sys.stderr.write("Bed attempting to read a line with less than three entries. The line looks like this: " + str(line))
        else:
            key = line[0]
            vals = [int(line[1]), int(line[2])]
            self.add_entry(key, vals)

    def read_file(self, reader):
        for line in reader:
            if line[0].startswith('#'):
                continue
            else:
                self.process_line(line)



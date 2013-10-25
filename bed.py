#!/usr/bin/env python

class Bed:


    def __init__(self, entries=None):
        if entries == None:
            self.entries = {}
        else:
            self.entries = entries

    def add_entry(self, key, val):
        self.entries[key] = val

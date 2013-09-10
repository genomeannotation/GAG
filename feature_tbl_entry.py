#!/usr/bin/env python

import copy

#TODO: It'd be cool if we could pickle/depickle
class FeatureTblEntry:

    #TODO: Add parameters to constructor
    def __init__(self):
        self.type = ''
        self.coords = []
        self.strand = '+'
        self.phase = 0 # Start codon of 1
        self.has_start = False
        self.has_stop = False
        self.annotations = []

    def set_type(self, newType):
        self.type = newType

    def add_coordinates(self, start, stop):
        if start >= stop:
            print("Warning: feature start >= stop; skipping...")
        else:
            self.coords.append([start, stop])

    def set_strand(self, strand):
        self.strand = strand

    def set_phase(self, phase):
        self.phase = phase

    def set_partial_info(self, has_start, has_stop):
        self.has_start = has_start
        self.has_stop = has_stop

    def add_annotation(self, key, value):
        self.annotations.append([key, value])

    # Returns the summed length of each coordinate pair
    def get_total_length(self):
        length = 0
        for coords in self.coords:
            length += coords[1]-coords[0]
        return length

    def write_to_string(self):
        entry = ''

        fixedCoords = copy.deepcopy(self.coords) # Make a deep copy of our list of coordinates so we don't mess with the original data
        fixedCoords.sort()
        if self.strand == '-':
            fixedCoords.reverse()
            [coords.reverse() for coords in fixedCoords]

        # Write the first pair of coords with the entry type and partial info
        if not self.has_start:
            entry += '<'
        entry += str(fixedCoords[0][0])+'\t'+str(fixedCoords[0][1])+'\t'+self.type+'\n'
        fixedCoords = fixedCoords[1:] # Cut out the first coordinates        

        if len(fixedCoords) > 0:
            # Write the middle pairs of coords
            for coords in fixedCoords[:-1]: # TODO: List comprehension?
                entry += str(coords[0])+'\t'+str(coords[1])+'\n'

            # Write the last pair of coords with partial info
            if self.has_stop:
                entry += str(fixedCoords[-1][0])+'\t'+str(fixedCoords[-1][1])+'\n'
            else:
                entry += str(fixedCoords[-1][0])+'\t'+'>'+str(fixedCoords[-1][1])+'\n'

        # Write the start codon annotation
        if self.phase != 0:
            entry += '\t\t\tstart_codon\t'+str(self.phase+1)+'\n'

        for annot in self.annotations:
            entry += '\t\t\t'+annot[0]+'\t'+annot[1]+'\n'

        return entry
        
        

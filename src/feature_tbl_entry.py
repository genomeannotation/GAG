#!/usr/bin/env python

import copy

class FeatureTblEntry:

    def __init__(self):
        self.type = ''
        self.name = ''
        self.seq_name = ''
        self.coords = []
        self.strand = '+'
        self.phase = 0 # Start codon of 1
        self.has_start = False
        self.has_stop = False
        self.annotations = []

    def set_type(self, newType):
        self.type = newType

    def set_name(self, name):
        self.name = name

    def set_seq_name(self, name):
        self.seq_name = name

    def add_coordinates(self, start, stop):
        if start >= stop:
            print("Warning: Entry "+self.name+": start >= stop")
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

    def add_annotations(self, annotations):
        for annotation in annotations:
            self.annotations.append(annotation)

    def get_annotation(self, key):
        for annotation in self.annotations:
            if annotation[0] == key:
                return annotation[1]
        return None

    def remove_annotation(self, key):
        to_remove = None
        for annotation in self.annotations:
            if annotation[0] == key:
                to_remove = annotation
        if to_remove:
            self.annotations.remove(to_remove)

    def is_hypothetical(self):
        for annotation in self.annotations:
            if annotation[0] == "product" and \
               annotation[1] == "hypothetical protein":
                return True
        return False

    def write_to_string(self):
        entry = ''

        if len(self.coords) == 0:
            print("Warning: Entry "+self.name+\
                  ": No coordinates; skipping entry")
            return entry

        # Make a deep copy of our list of coordinates 
        # so we don't mess with the original data
        fixedCoords = copy.deepcopy(self.coords) 
        fixedCoords.sort()
        if self.strand == '-':
            fixedCoords.reverse()
            [coords.reverse() for coords in fixedCoords]

        # Write the first pair of coords with the entry type and partial info
        if not self.has_start:
            entry += '<'
        if len(fixedCoords) == 1 and not self.has_stop: 
            # Special case: Only one set of coordinates and 
            # no stop codon - write the partial infos
            entry += str(fixedCoords[0][0])+'\t'+\
                     '>'+str(fixedCoords[0][1])+'\t'+self.type+'\n'
        else: 
            # Normal case - more than one pair of coordinates 
            # or we do have a stop codon
            entry += str(fixedCoords[0][0])+'\t'+\
                     str(fixedCoords[0][1])+'\t'+self.type+'\n'
        fixedCoords = fixedCoords[1:] # Cut out the first coordinates

        if len(fixedCoords) > 0:
            # Write the middle pairs of coords
            for coords in fixedCoords[:-1]: 
                entry += str(coords[0])+'\t'+str(coords[1])+'\n'

            # Write the last pair of coords with partial info
            if self.has_stop:
                entry += str(fixedCoords[-1][0])+\
                         '\t'+str(fixedCoords[-1][1])+'\n'
            else:
                entry += str(fixedCoords[-1][0])+'\t'+\
                         '>'+str(fixedCoords[-1][1])+'\n'

        # Write the start codon annotation
        if self.type == 'CDS':
            if self.phase != 0 and not self.has_start:
                entry += '\t\t\tcodon_start\t'+str(self.phase+1)+'\n'
            else:
                entry += '\t\t\tcodon_start\t1\n'

        for annot in self.annotations:
            if len(annot) == 2:
                entry += '\t\t\t'+annot[0]+'\t'+annot[1]+'\n'
            else:
                print("ERROR: Invalid annotation on "+self.name+": "+annot)

        return entry
        
        

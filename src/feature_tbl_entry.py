#!/usr/bin/env python

import copy

#TODO: It'd be cool if we could pickle/depickle
class FeatureTblEntry:

    #TODO: Add parameters to constructor
    def __init__(self):
        self.type = ''
        self.name = ''
        self.coords = []
        self.strand = '+'
        self.phase = 0 # Start codon of 1
        self.has_start = False
        self.has_stop = False
        self.annotations = []
        self.errors = [] # errors associated with this entry

    def set_type(self, newType):
        self.type = newType

    def set_name(self, name):
        self.name = name

    def add_error(self, error):
        self.errors.append(error)

    def has_error(self, error):
        for err in errors:
            if err == error:
                return True
        return False

    def add_coordinates(self, start, stop):
        if start >= stop:
            print("Warning: Entry "+self.name+": start >= stop; skipping coordinate")
        else:
            self.coords.append([start, stop])

    def has_same_coords(self, other):
        if len(self.coords) != len(other.coords):
            return False

        for i in range(len(self.coords)):
            if self.coords[i][0] != other.coords[i][0] or self.coords[i][1] != other.coords[i][1]:
                return False
        return True

    def is_short_intron(self):
        for coords in self.coords:
            if coords[1]-coords[0] <= 10: #NCBI says anything less than 10 nt is a short intron
                return True
        return False

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

    # Returns the summed length of each coordinate pair
    def get_total_length(self):
        length = 0
        for coords in self.coords:
            length += coords[1]-coords[0]
        return length

    def write_to_string(self):
        entry = ''

        if len(self.coords) == 0:
            print("Warning: Entry "+self.name+": No coordinates; skipping entry")
            return entry

        fixedCoords = copy.deepcopy(self.coords) # Make a deep copy of our list of coordinates so we don't mess with the original data
        fixedCoords.sort()
        if self.strand == '-':
            fixedCoords.reverse()
            [coords.reverse() for coords in fixedCoords]

        # Write the first pair of coords with the entry type and partial info
        if not self.has_start:
            entry += '<'
        if len(fixedCoords) == 1 and not self.has_stop: # Special case: Only one set of coordinates and no stop codon - write the partial infos
            entry += str(fixedCoords[0][0])+'\t'+'>'+str(fixedCoords[0][1])+'\t'+self.type+'\n'
        else: # Normal case - more than one pair of coordinates or we do have a stop codon
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
        
        

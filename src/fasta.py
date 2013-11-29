#!/usr/bin/env python

from bed import Bed

class Fasta:

    def __init__(self):
        self.entries = list()

    def read_string(self, data):
        seq_id = ''
        seq = ''
        for line in data.split('\n'):
            if line[0] == '>':
                if len(seq_id) > 0:
                    # Save the data
                    self.entries.append([seq_id, seq])
                
                seq_id = line[1:].strip().split()[0] # Get the next seq_id
                seq = ''
            else:
                seq += line.strip()
        # Add the last sequence
        self.entries.append([seq_id, seq])

    def write_string(self):
        s = str()
        i = 0
        for entry in self.entries:
            if i != 0:
                s += '\n'
            s += '>'+entry[0]+'\n'+entry[1]
            i += 1
        return s

    def get_seq(self, seq_id):
        for entry in self.entries:
            if entry[0] == seq_id:
                return entry[1]
        return None

    def get_subseq(self, seq_id, indices):
        seq = self.get_seq(seq_id)
        if seq and self.indices_not_out_of_range(indices, seq):
            start = indices[0] - 1
            end = indices[1]
            return seq[start:end]
        else:
            return None
        

    def read_file(self, fileName):
        with open(fileName, 'r') as f:
            seq_id = ''
            seq = ''
            for line in f:
                if line[0] == '>':
                    if len(seq_id) > 0:
                        # Save the data
                        self.entries.append([seq_id, seq])
                    seq_id = line[1:].strip().split()[0] # Get the next seq_id
                    seq = ''
                else:
                    seq += line.strip()
            # Add the last sequence
            self.entries.append([seq_id, seq])

    # i think this does the inverse of what we need?
    def trim_seq(self, seq_id, start, stop):
        for i, entry in enumerate(self.entries):
            if entry[0] == seq_id:
                if len(entry[1]) <= stop-(start-1):
                    self.entries.pop(i)
                else:
                    self.entries[i][1] = entry[1][:start-1] + entry[1][stop:]

    def apply_bed(self, bed):
        seqs_to_remove = []
        for i, seq in enumerate(self.entries):
            if bed.contains(seq[0]):
                coords = bed.get_coordinates(seq[0])
                if coords == [0, 0]:
                    seqs_to_remove.append(i)
                    seq[1] = coords
                else:
                    seq[1] = seq[1][coords[0]-1:coords[1]]
        # remove entries with coords = [0, 0]
        self.entries = [e for e in self.entries if e[1] != [0, 0]] 

    # 'indices' in question are one-based, not zero-based.
    def indices_not_out_of_range(self, indices, seq):
        if indices[0] > 0 and indices[1] <= len(seq):
            return True
        else:
            return False

#!/usr/bin/env python

class Fasta:

    def __init__(self):
        self.entries = list()

    def __str__(self):
        result = "Fasta containing "
        result += str(len(self.entries))
        result += " sequences\n"
        return result

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
        seq = ''

        for i in indices:
            subseq = self.get_seq(seq_id)
            if subseq and self.indices_not_out_of_range(i, subseq):
                start = i[0] - 1
                end = i[1]
                seq += subseq[start:end]
            else:
                return None
        return seq
        

    def read(self, io_buffer):
        seq_id = ''
        seq = ''
        for line in io_buffer:
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

    # 'indices' in question are one-based, not zero-based.
    def indices_not_out_of_range(self, indices, seq):
        if indices[0] > 0 and indices[1] <= len(seq):
            return True
        else:
            return False

    def subset_fasta(self, seqlist):
        self.entries = [e for e in self.entries if e[0] in seqlist]

    def remove_seq(self, seq_id):
        self.entries = [e for e in self.entries if e[0] != seq_id]

    # Given a position in the fasta, returns the number of Ns 
    # from that position forward 
    # (returns 0 if the base at that position is not N)
    def how_many_Ns_forward(self, seq_id, position):
        seq = self.get_seq(seq_id)
        index = position-1
        if seq[index] != 'N' and seq[index] != 'n':
            return 0
        else:
            count = 1
            index += 1
            for base in seq[index:]:
                if base != 'N' and base != 'n':
                    return count
                else:
                    count += 1
            return count

    # Given a position in the fasta, returns the number of Ns 
    # from that position backward 
    # (returns 0 if the base at that position is not N)
    def how_many_Ns_backward(self, seq_id, position):
        seq = self.get_seq(seq_id)
        index = position-1
        if seq[index] != 'N' and seq[index] != 'n':
            return 0
        else:
            count = 1
            index -= 1
            for base in seq[index::-1]:
                if base != 'N' and base != 'n':
                    return count
                else:
                    count += 1
            return count



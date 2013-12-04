#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import sys

bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

def valid_seq(seq):
    # Assumes seq is already lowercase
    for base in seq:
        if base not in ['a', 'c', 't', 'g', 'n']:
            return False
    if len(seq) < 3:
        return False
    else:
        return True

def valid_frame(frame):
    return frame in [1, 2, 3]

def valid_strand(strand):
    return strand in ['+', '-']

def verify_inputs(seq, frame, strand):
    if not valid_seq(seq):
        sys.stderr.write("Invalid seq passed to translate.py: " + seq + "\n")
        return False
    elif not valid_frame(frame):
        sys.stderr.write("Invalid frame passed to translate.py: " + str(frame) + "\n")
        sys.stderr.write("(Frame should be 1, 2 or 3)")
        return False
    elif not valid_strand(strand):
        sys.stderr.write("Invalid strand passed to translate.py: " + strand + "\n")
        return False
    else:
        return True

def reverse_complement(seq):
    bases = ['a', 'c', 'g', 't']
    complements = ['t', 'g', 'c', 'a']
    rev_comp_dict = dict(zip(bases, complements))
    return ''.join([rev_comp_dict.get(base) for base in reversed(seq)])

def translate(seq, frame, strand):
    # TODO handle N's?
    seq = seq.lower().replace('\n', '').replace(' ', '')

    if not verify_inputs(seq, frame, strand):
        return None

    # Adjust according to strand
    if strand == '-':
        seq = reverse_complement(seq)

    # Adjust according to frame
    seq = seq[frame-1:]

    # Now translate
    peptide = ''
    
    for i in xrange(0, len(seq), 3):
        codon = seq[i: i+3]
        amino_acid = codon_table.get(codon, '')
        peptide += amino_acid

    return peptide

#!/usr/bin/env python
# coding=utf-8
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4


BASES = ['t', 'c', 'a', 'g']
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AMINO_ACIDS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON_TABLE = dict(zip(CODONS, AMINO_ACIDS))


def valid_seq(seq):
    # Assumes seq is already lowercase
    for base in seq:
        if base not in ['a', 'c', 't', 'g', 'n']:
            return False
    if len(seq) < 3:
        return False
    else:
        return True


def valid_strand(strand):
    return strand in ['+', '-']


def has_start_codon(seq):
    return seq[:3].lower() == 'aug' or seq[:3].lower() == 'atg'


def has_stop_codon(seq):
    last3 = seq[-3:].lower()
    if last3 == 'tag':
        return True
    elif last3 == 'taa':
        return True
    elif last3 == 'tga':
        return True
    else:
        return False


def contains_internal_stop(seq, strand):
    translated = translate(seq, strand)
    return '*' in translated[:-1]


def reverse_complement(seq):
    complement = {'a': 't',
                  'A': 'T',
                  'c': 'g',
                  'G': 'C',
                  'g': 'c',
                  'C': 'G',
                  't': 'a',
                  'T': 'A',
                  'n': 'n',
                  'N': 'N'}
    return ''.join([complement.get(base, 'N') for base in seq[::-1]])


def translate(seq, strand):
    seq = seq.lower().replace('\n', '').replace(' ', '')

    # Verify strand
    if not valid_strand(strand):
        return ""

    # Adjust according to strand
    if strand == '-':
        seq = reverse_complement(seq)

    # Now translate
    peptide = ''

    for i in xrange(0, len(seq), 3):
        codon = seq[i: i + 3]
        if len(codon) != 3:
            amino_acid = ''
        elif 'N' in codon or 'n' in codon or codon not in CODON_TABLE.keys():
            amino_acid = 'X'
        else:
            amino_acid = CODON_TABLE.get(codon, '')
        peptide += amino_acid

    return peptide

#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

def translate(seq, frame, strand):
    seq = seq.lower().replace('\n', '').replace(' ', '')

    # Adjust according to frame and strand
    if strand == '-':
        seq = seq[len(seq)-1-(frame-1)::-1]
        lSeq = list(seq)
        for i in range(0, len(lSeq)):
            if lSeq[i] == 'a':
                lSeq[i] = 't'
            elif lSeq[i] == 't':
                lSeq[i] = 'a'
            elif lSeq[i] == 'g':
                lSeq[i] = 'c'
            elif lSeq[i] == 'c':
                lSeq[i] = 'g'
        seq = ''.join(lSeq)
    elif strand == '+':
        seq = seq[frame-1:]
    else:
        print('ERROR: Invalid strand\n')
        return str()

    print(seq+'\n')
    peptide = ''
    
    for i in xrange(0, len(seq), 3):
        codon = seq[i: i+3]
        amino_acid = codon_table.get(codon, '')
        peptide += amino_acid
                
    return peptide

seq = 'CGGCTCCAAAaAGTGCAACAAACACGAACAAAAACTTATTCAT'
print(translate(seq, 2, '-'))

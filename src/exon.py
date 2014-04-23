#!/usr/bin/env python

import copy
from src.gene_part import *

class Exon(GenePart):

    def __init__(self, **kwargs):
        kwargs['feature_type'] = 'exon'
        GenePart.__init__(self, **kwargs)

    def to_tbl(self, has_start, has_stop):
        indices = copy.deepcopy(self.indices)
        return write_tbl_entry(indices, self.strand, has_start, has_stop, False, self.annotations)


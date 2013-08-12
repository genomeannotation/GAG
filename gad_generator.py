#!/usr/bin/env python

from gff_reader import GffReader
from sqlite_wrapper import SqliteWrapper

gff = GffReader()
gff_db = gff.load('test_files/test.gff')
row = gff_db.getRowAsStr('gff', 'id', '1')
print(row)

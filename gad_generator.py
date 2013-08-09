#!/usr/bin/env python

from gff_reader import GffReader
from sqlite_wrapper import SqliteWrapper

gff = GffReader()
db = gff.load('test_files/test.gff')
sqlite = SqliteWrapper(db)
row = sqlite.getRowAsStr('gff', 'id', '1')
print(row)

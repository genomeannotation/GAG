#!/usr/bin/env python

from gff_reader import GffReader
from sqlite_wrapper import SqliteWrapper

gff = GffReader()
gff_db = gff.read_into_db('test_files/test.gff', ':memory:')

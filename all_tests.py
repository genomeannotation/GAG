#!/usr/bin/env python
# coding=utf-8

# import all the lovely files
import unittest
from test import fasta_reader_tests, gene_part_tests, xrna_tests, gene_tests, translator_tests, gff_reader_tests,\
    sequence_tests, filter_manager_tests, filters_tests, stats_manager_tests, seq_helper_tests, cds_tests, exon_tests

# get suites from test modules
suite1 = fasta_reader_tests.suite()
suite2 = gene_part_tests.suite()
suite3 = xrna_tests.suite()
suite4 = gene_tests.suite()
suite9 = translator_tests.suite()
suite10 = gff_reader_tests.suite()
suite11 = sequence_tests.suite()
suite12 = filter_manager_tests.suite()
suite13 = filters_tests.suite()
suite14 = stats_manager_tests.suite()
suite15 = seq_helper_tests.suite()
suite16 = cds_tests.suite()
suite17 = exon_tests.suite()

# collect suites in a TestSuite object
suite = unittest.TestSuite()
suite.addTest(suite1)
suite.addTest(suite2)
suite.addTest(suite3)
suite.addTest(suite4)
suite.addTest(suite9)
suite.addTest(suite10)
suite.addTest(suite11)
suite.addTest(suite12)
suite.addTest(suite13)
suite.addTest(suite14)
suite.addTest(suite15)
suite.addTest(suite16)
suite.addTest(suite17)

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)

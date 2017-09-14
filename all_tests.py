#!/usr/bin/env python
# coding=utf-8

# import all the lovely files
import unittest
from test import fasta_reader_tests, gene_part_tests, xrna_tests, gene_tests, translator_tests, gff_reader_tests,\
    sequence_tests, filter_manager_tests, filters_tests, stats_manager_tests, seq_helper_tests, cds_tests, exon_tests

# collect suites in a TestSuite object
suite = unittest.TestSuite()
suite.addTest(fasta_reader_tests.suite())
suite.addTest(gene_part_tests.suite())
suite.addTest(xrna_tests.suite())
suite.addTest(gene_tests.suite())
suite.addTest(translator_tests.suite())
suite.addTest(gff_reader_tests.suite())
suite.addTest(sequence_tests.suite())
suite.addTest(filter_manager_tests.suite())
suite.addTest(filters_tests.suite())
suite.addTest(stats_manager_tests.suite())
suite.addTest(seq_helper_tests.suite())
suite.addTest(cds_tests.suite())
suite.addTest(exon_tests.suite())

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)

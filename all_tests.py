#!/usr/bin/env python

# import all the lovely files
import unittest
import test.fasta_tests
import test.gene_part_tests
import test.mrna_tests
import test.gene_tests
import test.gff_tests
import test.assorted_tests
import test.genome_tests
import test.console_controller_tests
import test.annotator_tests
import test.translate_tests
import test.gff_reader_tests

# get suites from test modules
suite1 = test.fasta_tests.suite()
suite2 = test.gene_part_tests.suite()
suite3 = test.mrna_tests.suite()
suite4 = test.gene_tests.suite()
suite5 = test.gff_tests.suite()
suite6 = test.assorted_tests.suite()
suite7 = test.genome_tests.suite()
suite8 = test.console_controller_tests.suite()
suite9 = test.annotator_tests.suite()
suite10 = test.translate_tests.suite()
suite11 = test.gff_reader_tests.suite()

# collect suites in a TestSuite object
suite = unittest.TestSuite()
suite.addTest(suite1)
suite.addTest(suite2)
suite.addTest(suite3)
suite.addTest(suite4)
suite.addTest(suite5)
suite.addTest(suite6)
suite.addTest(suite7)
suite.addTest(suite8)
suite.addTest(suite9)
suite.addTest(suite10)
suite.addTest(suite11)

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)

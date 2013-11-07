#!/usr/bin/env python

# import all the lovely files
import unittest
import test.fasta_tests
import test.gene_part_tests
import test.mrna_tests
import test.gene_tests
import test.bed_tests
import test.gff_tests

# get suites from test modules
suite1 = test.fasta_tests.suite()
suite2 = test.gene_part_tests.suite()
suite3 = test.mrna_tests.suite()
suite4 = test.gene_tests.suite()
suite5 = test.bed_tests.suite()
suite6 = test.gff_tests.suite()

# collect suites in a TestSuite object
suite = unittest.TestSuite()
suite.addTest(suite1)
suite.addTest(suite2)
suite.addTest(suite3)
suite.addTest(suite4)
suite.addTest(suite5)
suite.addTest(suite6)

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)

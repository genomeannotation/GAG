#!/usr/bin/env python

import unittest
import test.fasta_tests
import test.feature_class_tests
import test.bed_tests
import test.gff_tests

# get suites from test modules
suite1 = test.feature_class_tests.suite()
suite2 = test.bed_tests.suite()
suite3 = test.fasta_tests.suite()
suite4 = test.gff_tests.suite()

# collect suites in a TestSuite object
suite = unittest.TestSuite()
suite.addTest(suite1)
suite.addTest(suite2)
suite.addTest(suite3)
suite.addTest(suite4)

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)

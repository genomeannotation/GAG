#!/usr/bin/env python

import unittest
import fasta_tests
import feature_class_tests
import bed_tests
import gff_tests

# get suites from test modules
suite1 = feature_class_tests.suite()
suite2 = bed_tests.suite()
suite3 = fasta_tests.suite()
suite4 = gff_tests.suite()

# collect suites in a TestSuite object
suite = unittest.TestSuite()
suite.addTest(suite1)
suite.addTest(suite2)
suite.addTest(suite3)
suite.addTest(suite4)

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)

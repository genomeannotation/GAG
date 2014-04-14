#!/usr/bin/env python

import unittest
from src.seq_fixer import SeqFixer

class TestSeqFixer(unittest.TestCase):

    def setUp(self):
        self.fixer = SeqFixer()

    def test_initialize(self):
        self.assertFalse(self.fixer.terminal_ns)
        self.assertFalse(self.fixer.internal_stops)
        self.assertFalse(self.fixer.start_stop_codons)


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSeqFixer))
    return suite

if __name__ == '__main__':
    unittest.main()

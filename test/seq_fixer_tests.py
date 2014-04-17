#!/usr/bin/env python

import unittest
from src.seq_fixer import SeqFixer

class TestSeqFixer(unittest.TestCase):

    def setUp(self):
        self.fixer = SeqFixer()

    def test_initialize(self):
        self.assertFalse(self.fixer.terminal_ns)
        self.assertFalse(self.fixer.start_stop_codons)

    def test_fix_terminal_ns(self):
        self.assertFalse(self.fixer.terminal_ns)
        self.fixer.fix_terminal_ns()
        self.assertTrue(self.fixer.terminal_ns)

    def test_fix_start_stop_codons(self):
        self.assertFalse(self.fixer.start_stop_codons)
        self.fixer.fix_start_stop_codons()
        self.assertTrue(self.fixer.start_stop_codons)

    def test_dirty(self):
        self.assertFalse(self.fixer.dirty)


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSeqFixer))
    return suite

if __name__ == '__main__':
    unittest.main()

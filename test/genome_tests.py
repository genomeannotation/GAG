#!/usr/bin/env python

import unittest
from src.fasta import Fasta
from src.genome import Genome

class TestGenome(unittest.TestCase):

    def setUp(self):
        self.genome = Genome()

    def test_constructor(self):
        self.assertEqual('Genome', self.genome.__class__.__name__)

    def test_add_template_file(self):
        self.genome.add_template_file("sample_files/sample.sbt")
        self.assertTrue(self.genome.template_file)
        self.assertEquals("sample_files/sample.sbt", self.genome.template_file)



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGenome))
    return suite

if __name__ == '__main__':
    unittest.main()

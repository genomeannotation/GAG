#!/usr/bin/env python

import unittest
import io
from src.fasta import Fasta
from mock import Mock

class TestFasta(unittest.TestCase):

    def setUp(self):
        self.fasta0 = Fasta()
        self.fasta1 = Fasta()
        self.fasta1.entries = [['seq1', 'ATGCCGTA'], ['seq2', 'AGGTCC'], ['seq3', 'GGGGGG']]

    def test_constructor(self):
        self.assertEquals('Fasta', self.fasta0.__class__.__name__)
        self.assertTrue(isinstance(self.fasta0.entries, list))

    def test_string(self):
        expected = "Fasta containing 3 sequences\n"
        self.assertEquals(expected, str(self.fasta1))

    def test_trim_seq(self):
        self.fasta1.trim_seq('seq2', 2, 3)
        expected = '>seq1\nATGCCGTA\n>seq2\nATCC\n>seq3\nGGGGGG'
        self.assertEqual(expected, self.fasta1.write_string())
        self.fasta1.trim_seq('seq2', 1, 4)
        expected = '>seq1\nATGCCGTA\n>seq3\nGGGGGG'
        self.assertEqual(expected, self.fasta1.write_string())

    def test_get_seq(self):
        expected = 'AGGTCC'
        self.assertEqual(expected, self.fasta1.get_seq('seq2'))
        self.assertFalse(self.fasta1.get_seq('bad_seqid'))

    def test_get_subseq(self):
        expected = 'GGTC'
        self.assertEqual(expected, self.fasta1.get_subseq('seq2', [[2, 5]]))
        self.assertFalse(self.fasta1.get_subseq('bad_seqid', [[1, 10]]))
        self.assertFalse(self.fasta1.get_subseq('seq1', [[1, 10]])) #out of range

    def test_read(self):
        no_line_breaks = io.BytesIO('>seq_1\nGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA\n' +
                                    '>seq_2\nNNNNNNNNGATTACAGATTACAGATTACANNNNNNNNNNN')
        line_breaks = io.BytesIO('>seq_1\nGATTACAGATTACAGATTACAGATTACA\nGATTACAGATTACAGATTACAGATTACA\n' +
                                 '>seq_2\nNNNNNNNNGATTACAGATTACAGATTAC\nANNNNNNNNNNN')

        self.fasta0.read(no_line_breaks)
        self.assertEquals(2, len(self.fasta0.entries))
        self.assertEquals('seq_1', self.fasta0.entries[0][0])
        self.assertEquals('GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA', self.fasta0.entries[0][1])
        self.assertEquals('seq_2', self.fasta0.entries[1][0])
        self.assertEquals('NNNNNNNNGATTACAGATTACAGATTACANNNNNNNNNNN', self.fasta0.entries[1][1])
        self.fasta0.read(line_breaks)
        self.assertEquals(4, len(self.fasta0.entries))
        self.assertEquals(['seq_1', 'GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA'], self.fasta0.entries[2])
        self.assertEquals(['seq_2', 'NNNNNNNNGATTACAGATTACAGATTACANNNNNNNNNNN'], self.fasta0.entries[3])

    def test_indices_not_out_of_range(self):
        test_seq = 'GATTACA'
        good_indices = [1, 7]
        bad_indices = [2, 8]
        self.assertTrue(self.fasta1.indices_not_out_of_range(good_indices, test_seq))

    def test_subset_fasta(self):
        self.assertEquals(3, len(self.fasta1.entries))
        seqs_to_keep = ['seq1', 'seq3']
        self.fasta1.subset_fasta(seqs_to_keep)
        self.assertEquals(2, len(self.fasta1.entries))
        self.assertEquals('seq1', self.fasta1.entries[0][0])
        self.assertEquals('seq3', self.fasta1.entries[1][0])

    def test_remove_seq(self):
        self.assertEqual(3, len(self.fasta1.entries))
        self.fasta1.remove_seq('seq1')
        self.assertEqual(2, len(self.fasta1.entries))
        self.assertEqual('seq2', self.fasta1.entries[0][0])

    def setup_bad_fasta(self):
        self.badfasta = Fasta()
        self.badfasta.entries = [['seq1', 'NNnNNGATTACA'], ['seq2', 'GATTACAnNn']]

    def test_how_many_Ns_forward(self):
        self.setup_bad_fasta()
        self.assertEqual(5, self.badfasta.how_many_Ns_forward('seq1', 1))

    def test_how_many_Ns_forward_returns_zero_if_no_Ns(self):
        # (seq1, base 3 is a 'G')
        self.assertEqual(0, self.fasta1.how_many_Ns_forward('seq1', 3))

    def test_how_many_Ns_backward(self):
        self.setup_bad_fasta()
        self.assertEqual(3, self.badfasta.how_many_Ns_backward('seq2', 10))

    def test_how_many_Ns_backward_returns_zero_if_no_Ns(self):
        self.assertEqual(0, self.fasta1.how_many_Ns_backward('seq1', 3))

    def test_contains_seq_id(self):
        self.assertTrue(self.fasta1.contains_seq_id("seq1"))

        


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFasta))
    return suite

if __name__ == '__main__':
    unittest.main()

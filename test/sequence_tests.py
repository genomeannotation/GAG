#!/usr/bin/env python

import unittest
from mock import Mock
from src.sequence import Sequence

class TestSequence(unittest.TestCase):

    def setUp(self):
        self.seq1 = Sequence("seq1", "GATTACA")

    def add_mock_gene(self):
        mockgene = Mock()
        mockgene.identifier = "foo_gene"
        mockgene.death_flagged = False
        mockgene.get_valid_mrnas = Mock(return_value=[])
        self.seq1.add_gene(mockgene)
        
    def add_mock_gene_with_1_mrna(self, name):
        mockgene = Mock()
        mockgene.identifier = name
        mockgene.death_flagged = False
        mockgene.mrnas = [Mock()]
        mockgene.mrnas[0].identifier = name+"-RA"
        mockgene.mrnas[0].cds = Mock()
        mockgene.mrnas[0].cds.identifier = [name+"-RA:CDS"]
        mockgene.mrnas[0].cds.length = Mock(return_value=5)
        mockgene.mrnas[0].exon = Mock()
        mockgene.mrnas[0].length = Mock(return_value=2)
        mockgene.get_valid_mrnas = Mock(return_value=mockgene.mrnas)
        mockgene.length = Mock(return_value=20)
        mockgene.get_partial_info.return_value = {"complete": 1, "start_no_stop": 0, "stop_no_start": 1, "no_stop_no_start": 1}
        mockgene.get_num_exons.return_value = 5
        mockgene.get_num_introns.return_value = 4
        mockgene.get_longest_exon.return_value = 20
        mockgene.get_longest_intron.return_value = 20
        mockgene.get_shortest_exon.return_value = 8
        mockgene.get_shortest_intron.return_value = 8
        mockgene.get_total_exon_length.return_value = 15
        mockgene.get_total_intron_length.return_value = 15
        self.seq1.add_gene(mockgene)
        
    def add_mock_gene_with_2_mrnas(self, name):
        mockgene = Mock()
        mockgene.identifier = name
        mockgene.death_flagged = False
        mockgene.mrnas = [Mock(), Mock()]
        mockgene.mrnas[0].identifier = name+"-RA"
        mockgene.mrnas[0].cds = None
        mockgene.mrnas[0].exon = None
        mockgene.mrnas[0].length = Mock(return_value=5)
        mockgene.mrnas[1].identifier = name+"-RB"
        mockgene.mrnas[1].cds = Mock()
        mockgene.mrnas[1].cds.identifier = [name+"-RB:CDS"]
        mockgene.mrnas[1].cds.length = Mock(return_value=3)
        mockgene.mrnas[1].exon = Mock()
        mockgene.mrnas[1].length = Mock(return_value=2)
        mockgene.get_valid_mrnas = Mock(return_value=mockgene.mrnas)
        mockgene.length = Mock(return_value=10)
        mockgene.get_partial_info.return_value = {"complete": 0, "start_no_stop": 1, "stop_no_start": 1, "no_stop_no_start": 1}
        mockgene.get_num_exons.return_value = 4
        mockgene.get_num_introns.return_value = 3
        mockgene.get_longest_exon.return_value = 10
        mockgene.get_longest_intron.return_value = 10
        mockgene.get_shortest_exon.return_value = 5
        mockgene.get_shortest_intron.return_value = 5
        mockgene.get_total_exon_length.return_value = 25
        mockgene.get_total_intron_length.return_value = 25
        self.seq1.add_gene(mockgene)

    def test_string(self):
        expected = "Sequence seq1 of length 7 containing 0 genes\n"
        self.assertEquals(expected, str(self.seq1))

    def test_how_many_Ns_forward(self):
        badseq = Sequence('seq1', 'NNnNNGATTACA')
        self.assertEqual(5, badseq.how_many_Ns_forward(1))

    def test_how_many_Ns_forward_returns_zero_if_no_Ns(self):
        badseq = Sequence('seq2', 'GATTACA')
        self.assertEqual(0, badseq.how_many_Ns_forward(3))

    def test_how_many_Ns_backward(self):
        badseq = Sequence('seq3', 'gattaNnN')
        self.assertEqual(3, badseq.how_many_Ns_backward(8))

    def test_how_many_Ns_backward_returns_zero_if_no_Ns(self):
        self.assertEqual(0, self.seq1.how_many_Ns_backward(3))

    def test_add_gene(self):
        self.add_mock_gene()
        self.assertEqual(1, len(self.seq1.genes))

    def test_remove_gene(self):
        self.add_mock_gene()
        self.assertEqual(1, len(self.seq1.genes))
        self.seq1.remove_gene("foo_gene")
        self.assertEqual(0, len(self.seq1.genes))

    def test_remove_gene_fails_if_no_match(self):
        self.add_mock_gene()
        self.assertEqual(1, len(self.seq1.genes))
        self.seq1.remove_gene("bar_gene")
        self.assertEqual(1, len(self.seq1.genes))

    def test_trim_region(self):
        self.assertEquals("GATTACA", self.seq1.bases)
        self.seq1.trim_region(1, 4)
        self.assertEquals("ACA", self.seq1.bases)

    def test_trim_region_trims_gene(self):
        # TODO don't remove gene, trim its indices and verify?
        self.add_mock_gene()
        self.seq1.genes[0].indices = [3, 7]
        self.assertEquals(1, len(self.seq1.genes))
        self.seq1.trim_region(1, 4)
        self.assertEquals(0, len(self.seq1.genes))

    def test_get_subseq(self):
        self.assertEquals("ATTA", self.seq1.get_subseq(2, 5))

    def test_get_cds_partial_info(self):
        self.add_mock_gene_with_1_mrna("foo_gene1")
        self.add_mock_gene_with_2_mrnas("foo_gene2")
        partial_info = self.seq1.get_cds_partial_info()
        self.assertEquals(1, partial_info["CDS: complete"])

    def test_to_tbl(self):
        self.add_mock_gene()
        self.seq1.genes[0].to_tbl.return_value = "mockgene to tbl"
        tbl = self.seq1.to_tbl()
        expected = ">Feature seq1\n"
        expected += "1\t7\tREFERENCE\n"
        expected += "\t\t\tPBARC\t12345\n"
        expected += "mockgene to tbl"
        self.assertEquals(tbl, expected)

    def test_stats(self):
        self.add_mock_gene_with_1_mrna("foo_gene1")
        self.add_mock_gene_with_2_mrnas("foo_gene2")
        stats = self.seq1.stats()
        self.assertEquals(stats["Total sequence length"], 7)
        self.assertEquals(stats["Number of genes"], 2)
        self.assertEquals(stats["Number of mRNAs"], 3)
        self.assertEquals(stats["Number of exons"], 9)
        self.assertEquals(stats["Number of introns"], 7)
        self.assertEquals(stats["Number of CDS"], 2)
        self.assertEquals(stats["CDS: complete"], 1)
        self.assertEquals(stats["CDS: start, no stop"], 1)
        self.assertEquals(stats["CDS: stop, no start"], 2)
        self.assertEquals(stats["CDS: no stop, no start"], 2)
        self.assertEquals(stats["Longest gene"], 20)
        self.assertEquals(stats["Longest mRNA"], 5)
        self.assertEquals(stats["Longest exon"], 20)
        self.assertEquals(stats["Longest intron"], 20)
        self.assertEquals(stats["Longest CDS"], 5)
        self.assertEquals(stats["Shortest gene"], 10)
        self.assertEquals(stats["Shortest mRNA"], 2)
        self.assertEquals(stats["Shortest exon"], 5)
        self.assertEquals(stats["Shortest intron"], 5)
        self.assertEquals(stats["Shortest CDS"], 3)
        self.assertEquals(stats["Total gene length"], 30)
        self.assertEquals(stats["Total mRNA length"], 9)
        self.assertEquals(stats["Total exon length"], 40)
        self.assertEquals(stats["Total intron length"], 40)
        self.assertEquals(stats["Total CDS length"], 8)
        
        self.seq1.genes[0].death_flagged = True
        stats = self.seq1.stats()
        self.assertEquals(stats["Number of genes"], 1)



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSequence))
    return suite

if __name__ == '__main__':
    unittest.main()

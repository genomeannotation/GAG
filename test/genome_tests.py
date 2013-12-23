#!/usr/bin/env python

import unittest
from src.fasta import Fasta
from src.genome import Genome
from mock import Mock

class TestGenome(unittest.TestCase):

    def setUp(self):
        self.genome = Genome()

    def test_constructor(self):
        self.assertEqual('Genome', self.genome.__class__.__name__)

    def test_verify_file(self):
        self.assertTrue(self.genome.verify_file("sample_files/sample.sbt"))
        self.assertFalse(self.genome.verify_file("no_file_here.foo"))

    def test_remove_all_gene_segments(self):
        self.assertFalse(self.genome.gff)
        gff = Mock()
        self.genome.gff = gff
        self.assertTrue(self.genome.gff)
        self.genome.remove_all_gene_segments("BDOR_foo")
        gff.remove_all_gene_segments.assert_called_with("BDOR_foo")

    def test_verify_start_codon(self):
        # this is a pretty mocky test. basically verifies that genome
        # gets cds indices from the mrna, gets corresponding subseq from fasta
        # and calls mrna's add_start_codon method
        mrna = Mock()
        mrna.get_cds_indices.return_value = [[100, 109], [150, 200]]
        fasta = Mock()
        fasta.get_subseq.return_value = 'auggattaca'
        self.genome.fasta = fasta
        seq_name = 'seq_1'
        self.genome.verify_start_codon(mrna, seq_name)
        mrna.get_cds_indices.assert_called_with()
        fasta.get_subseq.assert_called_with('seq_1', [100, 109])
        mrna.add_start_codon.assert_called_with(100)

    def test_verify_stop_codon(self):
        mrna = Mock()
        mrna.get_cds_indices.return_value = [[100, 109], [150, 200]]
        fasta = Mock()
        fasta.get_subseq.return_value = 'gattacatag'
        self.genome.fasta = fasta
        seq_name = 'seq_1'
        self.genome.verify_stop_codon(mrna, seq_name)
        mrna.get_cds_indices.assert_called_with()
        fasta.get_subseq.assert_called_with('seq_1', [150, 200])
        mrna.add_stop_codon.assert_called_with(200)

    def test_obliterate_genes_related_to_mrnas(self):
        gff = Mock()
        self.genome.gff = gff
        self.genome.obliterate_genes_related_to_mrnas(['foo', 'bar'])
        gff.obliterate_genes_related_to_mrnas.assert_called_with(['foo', 'bar'])

    def setup_mocks_for_rename_maker_mrnas(self):
        self.maker_mrna = Mock()
        self.maker_mrna.name = 'maker-guy'
        self.maker_mrna.is_maker_mrna.return_value = True
        self.bdor_mrna = Mock()
        self.bdor_mrna.name = 'bdor-guy'
        self.bdor_mrna.is_maker_mrna.return_value = False
        self.gene1 = Mock()
        self.gene1.mrnas = [self.maker_mrna, self.bdor_mrna]
        self.other_mrna = Mock()
        self.other_mrna.is_maker_mrna.return_value = False
        self.other_maker = Mock()
        self.other_maker.is_maker_mrna.return_value = True
        self.gene2 = Mock()
        self.gene2.mrnas = [self.other_mrna, self.other_maker]
        self.gff = Mock()
        self.gff.genes = [self.gene1, self.gene2]
        self.genome.gff = self.gff
        self.genome.annot = Mock()

    def test_rename_maker_mrnas(self):
        self.setup_mocks_for_rename_maker_mrnas()
        self.assertEquals('maker-guy', self.maker_mrna.name)
        self.genome.rename_maker_mrnas()
        self.assertEquals('BDOR_1000000', self.maker_mrna.name)
        self.assertEquals('bdor-guy', self.bdor_mrna.name)
        self.assertEquals('BDOR_1000001', self.other_maker.name)
        self.bdor_mrna.is_maker_mrna.assert_called_with()
        self.other_mrna.is_maker_mrna.assert_called_with()


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGenome))
    return suite

if __name__ == '__main__':
    unittest.main()

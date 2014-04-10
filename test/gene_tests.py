#!/usr/bin/env python

import unittest
from mock import Mock, PropertyMock
from src.gene_part import GenePart, CDS, Exon
from src.mrna import MRNA
from src.gene import Gene

class TestGene(unittest.TestCase):

    def setUp(self):
        self.test_gene0 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[3734, 7436], strand='+', identifier=1)
        self.test_gene1 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[3734, 7436], strand='+', identifier=1)
        self.fake_mrna1 = Mock()
        self.fake_mrna1.death_flagged = False
        self.fake_mrna2 = Mock()
        self.fake_mrna2.death_flagged = False
        self.test_gene1.mrnas.append(self.fake_mrna1)
        self.test_gene1.mrnas.append(self.fake_mrna2)

    def test_constructor(self):
        self.assertEqual('Gene', self.test_gene0.__class__.__name__)

    def test_length(self):
        self.assertEqual(3703, self.test_gene0.length())

    def test_length_of_shortest_cds_segment(self):
        self.fake_mrna1.length_of_shortest_cds_segment.return_value = 358
        self.fake_mrna2.length_of_shortest_cds_segment.return_value = 241
        self.assertEquals(241, self.test_gene1.length_of_shortest_cds_segment())
        self.fake_mrna1.length_of_shortest_cds_segment.assert_called_with()
        self.fake_mrna2.length_of_shortest_cds_segment.assert_called_with()

    def test_get_longest_exon(self):
        self.fake_mrna1.get_longest_exon.return_value = 10
        self.fake_mrna2.get_longest_exon.return_value = 20
        self.assertEquals(20, self.test_gene1.get_longest_exon())

    def test_get_shortest_exon(self):
        self.fake_mrna1.get_shortest_exon.return_value = 5
        self.fake_mrna2.get_shortest_exon.return_value = 8
        self.assertEquals(5, self.test_gene1.get_shortest_exon())

    def test_get_total_exon_length(self):
        self.fake_mrna1.get_total_exon_length.return_value = 15
        self.fake_mrna2.get_total_exon_length.return_value = 25
        self.assertEquals(40, self.test_gene1.get_total_exon_length())

    def test_get_num_exons(self):
        self.fake_mrna1.get_num_exons.return_value = 5
        self.fake_mrna2.get_num_exons.return_value = 4
        self.assertEquals(9, self.test_gene1.get_num_exons())

    def test_get_longest_intron(self):
        self.fake_mrna1.get_longest_intron.return_value = 8
        self.fake_mrna2.get_longest_intron.return_value = 10
        self.assertEquals(10, self.test_gene1.get_longest_intron())

    def test_get_shortest_intron(self):
        self.fake_mrna1.get_shortest_intron.return_value = 5
        self.fake_mrna2.get_shortest_intron.return_value = 8
        self.assertEquals(5, self.test_gene1.get_shortest_intron())

    def test_get_total_intron_length(self):
        self.fake_mrna1.get_total_intron_length.return_value = 15
        self.fake_mrna2.get_total_intron_length.return_value = 25
        self.assertEquals(40, self.test_gene1.get_total_intron_length())

    def test_get_num_introns(self):
        self.fake_mrna1.get_num_introns.return_value = 3
        self.fake_mrna2.get_num_introns.return_value = 2
        self.assertEquals(5, self.test_gene1.get_num_introns())

    def test_trim_region_does_nothing_when_region_is_after_gene(self):
        self.test_gene0 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[20, 40], strand='+', identifier=1)
        self.test_gene0.trim_region(45, 50)
        self.assertEquals([20, 40], self.test_gene0.indices)

    def test_trim_region_adjusts_indices_when_region_is_before_gene(self):
        self.test_gene0 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[20, 40], strand='+', identifier=1)
        self.test_gene0.trim_region(11, 15)
        self.assertEquals([15, 35], self.test_gene0.indices)

    def test_trim_region_handles_child_mrnas(self):
        self.test_gene0 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[20, 40], strand='+', identifier=1)
        self.test_gene0.mrnas = [Mock()]
        self.test_gene0.trim_region(11, 15)
        self.test_gene0.mrnas[0].adjust_indices.assert_called_with(-5)

    def test_trim_region_adjusts_indices_correctly_when_region_overlaps_gene(self):
        self.test_gene0 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[20, 40], strand='+', identifier=1)
        self.test_gene0.trim_region(16, 25)
        self.assertEquals([16, 30], self.test_gene0.indices)

    def test_trim_region_adjusts_indices_correctly_when_region_overlaps_gene_end(self):
        self.test_gene0 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[20, 40], strand='+', identifier=1)
        self.test_gene0.trim_region(35, 45)
        self.assertEquals([20, 34], self.test_gene0.indices)

    def test_trim_region_adjusts_indices_correctly_when_region_is_contained_inside_gene(self):
        self.test_gene0 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[20, 40], strand='+', identifier=1)
        self.test_gene0.trim_region(25, 35)
        self.assertEquals([20, 29], self.test_gene0.indices)

    def test_get_partial_info(self):
        self.fake_mrna1.has_stop.return_value = True
        self.fake_mrna1.has_start.return_value = True
        self.fake_mrna2.has_stop.return_value = False
        self.fake_mrna2.has_start.return_value = True
        results = self.test_gene1.get_partial_info()
        self.assertEquals(1, results["complete"])

    def test_adjust_indices(self):
        self.test_gene1.adjust_indices(16)
        self.fake_mrna1.adjust_indices.assert_called_with(16, 1)
        self.assertEquals(3750, self.test_gene1.indices[0])
        # adjust them back
        self.test_gene1.adjust_indices(-16)
        self.fake_mrna1.adjust_indices.assert_called_with(-16, 1)
        self.assertEquals(3734, self.test_gene1.indices[0])

    def test_invalidate_region_calls_cds_and_exon(self):
        self.test_gene1.mrnas[0].cds = Mock()
        self.test_gene1.mrnas[0].exon = Mock()
        self.test_gene1.mrnas[1].cds = Mock()
        self.test_gene1.mrnas[1].exon = Mock()
        self.test_gene1.invalidate_region(50, 100)
        self.test_gene1.mrnas[0].cds.invalidate_region.assert_called_with(50, 100)
        self.test_gene1.mrnas[0].exon.invalidate_region.assert_called_with(50, 100)
        self.test_gene1.mrnas[1].cds.invalidate_region.assert_called_with(50, 100)
        self.test_gene1.mrnas[1].exon.invalidate_region.assert_called_with(50, 100)

    def test_invalidate_region_beginning(self):
        self.test_gene1.invalidate_region(3730, 3737)
        self.assertEquals(self.test_gene1.indices, [3738, 7436])

    def test_invalidate_region_end(self):
        self.test_gene1.invalidate_region(7431, 7440)
        self.assertEquals(self.test_gene1.indices, [3734, 7430])

    def test_to_gff(self):
        self.fake_mrna1.to_gff.return_value = "fake mrna1 to gff here:)\n"
        self.fake_mrna2.to_gff.return_value = "fake mrna2 to gff here:)\n"
        expected = "sctg_0080_0020\tmaker\tgene\t3734\t7436\t.\t+\t."
        expected += "\tID=1;foo=dog\n"
        expected += "fake mrna1 to gff here:)\n"
        expected += "fake mrna2 to gff here:)\n"
        self.test_gene1.add_annotation('foo', 'dog')
        self.assertEquals(expected, self.test_gene1.to_gff())

    def test_str(self):
        expected = "Gene (ID=1, seq_name=sctg_0080_0020) containing 2 mrnas"
        self.assertEquals(expected, str(self.test_gene1))

    def test_create_starts_and_stops(self):
        mrna1 = Mock()
        mrna2 = Mock()
        self.test_gene0.mrnas = [mrna1, mrna2]
        seq_object = Mock()
        self.test_gene0.create_starts_and_stops(seq_object)
        mrna1.create_start_and_stop_if_necessary.assert_called_with(seq_object, '+')
        mrna2.create_start_and_stop_if_necessary.assert_called_with(seq_object, '+')

    def test_remove_first_cds_segment_if_shorter_than(self):
        self.test_gene1.remove_first_cds_segment_if_shorter_than(4)
        self.fake_mrna1.remove_first_cds_segment_if_shorter_than.assert_called_with(4)
        self.fake_mrna2.remove_first_cds_segment_if_shorter_than.assert_called_with(4)
        
    def test_to_tbl_positive(self):
        gene = Gene(seq_name="seq1", source="maker", indices=[1, 50], strand="+", identifier="foo_gene_1")
        self.assertFalse(gene.annotations)
        gene.add_annotation('foo', 'dog')
        mrna1 = Mock()
        mrna1.to_tbl.return_value = "mrna1_to_tbl...\n"
        mrna2 = Mock()
        mrna2.to_tbl.return_value = "mrna2_to_tbl...\n"
        gene.mrnas.append(mrna1)
        gene.mrnas.append(mrna2)
        expected = "1\t50\tgene\n\t\t\tlocus_tag\tfoo_gene_1\n\t\t\tfoo\tdog\nmrna1_to_tbl...\nmrna2_to_tbl...\n"
        self.assertEquals(gene.to_tbl(), expected)

    def test_gene_initialized_without_annotations(self):
        newgene = Gene(seq_name="seq1", source="maker", indices=[1, 50], strand="+", identifier="foo_gene_1")
        self.assertFalse(newgene.annotations)
        self.assertEquals(0, len(newgene.annotations))

    def test_gene_initialized_with_annotations(self):
        newgene = Gene(seq_name="seq1", source="maker", indices=[1, 50], strand="+", identifier="foo_gene_1",\
                annotations=[["bar", "cat"]])
        self.assertTrue(newgene.annotations)
        self.assertEquals(1, len(newgene.annotations))

    def test_to_tbl_negative(self):
        gene = Gene("seq1", "maker", [1, 50], "-", "foo_gene_1")
        mrna1 = Mock()
        mrna1.to_tbl.return_value = "mrna1_to_tbl...\n"
        mrna2 = Mock()
        mrna2.to_tbl.return_value = "mrna2_to_tbl...\n"
        gene.mrnas.append(mrna1)
        gene.mrnas.append(mrna2)
        expected = "50\t1\tgene\n\t\t\tlocus_tag\tfoo_gene_1\nmrna1_to_tbl...\nmrna2_to_tbl...\n"
        self.assertEquals(gene.to_tbl(), expected)

        


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGene))
    return suite

if __name__ == '__main__':
    unittest.main()

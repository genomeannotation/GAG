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
        self.fake_mrna2 = Mock()
        self.test_gene1.add_mrna(self.fake_mrna1)
        self.test_gene1.add_mrna(self.fake_mrna2)

    def test_constructor(self):
        self.assertEqual('Gene', self.test_gene0.__class__.__name__)

    def test_is_empty(self):
        self.assertTrue(self.test_gene0.is_empty())
        for mrna in self.test_gene1.mrnas:
            mrna.death_flagged = False
        self.assertFalse(self.test_gene1.is_empty())
        for mrna in self.test_gene1.mrnas:
            mrna.death_flagged = True
        self.assertTrue(self.test_gene1.is_empty())

    def test_length(self):
        self.assertEqual(3703, self.test_gene0.get_length())

    def test_add_mrna(self):
        self.assertEquals(0, len(self.test_gene0.mrnas))
        self.test_gene0.add_mrna(self.fake_mrna1)
        self.assertEquals(1, len(self.test_gene0.mrnas))
        self.test_gene0.add_mrna(self.fake_mrna2)
        self.assertEquals(2, len(self.test_gene0.mrnas))

    def test_contains_mrna_with_id(self):
        self.fake_mrna1.identifier = "BDOR_foo"
        self.assertTrue(self.test_gene1.contains_mrna_with_id("BDOR_foo"))
        self.assertFalse(self.test_gene1.contains_mrna_with_id("no_such_mrna_id"))

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

    def test_trim_end(self):
        self.test_gene1.trim_end(7400)
        self.assertEquals(3734, self.test_gene1.indices[0])
        self.assertEquals(7400, self.test_gene1.indices[1])
        # verify recursive call on child mrnas
        self.fake_mrna1.trim_end.assert_called_with(7400)
        self.fake_mrna2.trim_end.assert_called_with(7400)
        # if trim before beginning of gene, mark for removal
        # by setting indices = [0, 0]
        self.test_gene1.trim_end(50)
        self.assertEquals(0, self.test_gene1.indices[0])
        self.assertEquals(0, self.test_gene1.indices[1])

    def test_trim_begin(self):
        self.test_gene1.trim_begin(10)
        self.assertEquals(3725, self.test_gene1.indices[0])
        self.assertEquals(7427, self.test_gene1.indices[1])
        # should call adjust_indices on child mrnas ...
        self.fake_mrna1.adjust_indices.assert_called_with(-9, 1)
        self.fake_mrna2.adjust_indices.assert_called_with(-9, 1)
        
    def test_clean_up_indices(self):
        # if indices[0] < 1, set to 1 
        nice_gene = Gene(seq_name='fooseq', source='maker', indices=[-23, 127], strand='-', identifier='foo')
        mrna1 = Mock()
        mrna2 = Mock()
        nice_gene.add_mrna(mrna1)
        nice_gene.add_mrna(mrna2)
        nice_gene.clean_up_indices()
        self.assertEquals(1, nice_gene.indices[0])
        self.assertEquals(127, nice_gene.indices[1])
        # test recursive call on child mrnas
        mrna1.clean_up_indices.assert_called_with()
        mrna2.clean_up_indices.assert_called_with()

        # if indices [1] < 1, mark for removal by setting
        # indices = [0, 0]
        junk_gene = Gene(seq_name='barseq', source='maker', indices=[-400, -100], strand='-', identifier='bar')
        junk_gene.clean_up_indices()
        self.assertEquals(0, junk_gene.indices[0])
        self.assertEquals(0, junk_gene.indices[1])

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
        
    def test_remove_invalid_features(self):
        gene = Gene(seq_name="sctg_foo", source='maker', indices=[100, 200], strand='-', identifier='foo_gene')
        nice_mrna = Mock()
        nice_ind = PropertyMock(return_value = [100, 200])
        type(nice_mrna).indices = nice_ind
        junk_mrna = MRNA(identifier='junk', indices=[0, 0], parent_id='foo_gene')
        gene.add_mrna(junk_mrna)
        gene.add_mrna(nice_mrna) 
        self.assertEquals(2, len(gene.mrnas))
        # should remove junk_mrna, keeping nice_mrna
        gene.remove_invalid_features()
        self.assertEquals(1, len(gene.mrnas))
        # should call recursively on any mrnas not discarded
        nice_mrna.remove_invalid_features.assert_called_with()

    def test_trim(self):
        gene = Gene(seq_name='sctg_foo', source='maker', indices=[100, 200], strand='-', identifier='foo_gene')
        # set up mocks...
        nice_mrna = Mock()
        nice_ind = PropertyMock(return_value = [100, 200])
        type(nice_mrna).indices = nice_ind
        bad_mrna = Mock()
        bad_ind = PropertyMock(return_value = [0, 0])
        type(bad_mrna).indices = bad_ind
        # add mocks...
        gene.add_mrna(nice_mrna)
        gene.add_mrna(bad_mrna)
        # starting off with two mrnas
        self.assertEquals(2, len(gene.mrnas))
        # trim...
        bed_indices = [120, 180]
        gene.trim(bed_indices)
        # verify new indices for gene
        self.assertEquals(1, gene.indices[0])
        self.assertEquals(61, gene.indices[1])
        # verify appropriate calls were made to child mrnas
        nice_mrna.trim_end.assert_called_with(180)
        bad_mrna.trim_end.assert_called_with(180)
        nice_mrna.adjust_indices.assert_called_with(-119, 1)
        bad_mrna.adjust_indices.assert_called_with(-119, 1)
        nice_mrna.adjust_phase.assert_called_with()
        bad_mrna.adjust_phase.assert_called_with()
        nice_mrna.clean_up_indices.assert_called_with()
        bad_mrna.clean_up_indices.assert_called_with()
        nice_mrna.remove_invalid_features.assert_called_with()
        # verify bad_mrna was removed
        self.assertEquals(1, len(gene.mrnas))

        # if indices=[0, 0] then remove all features
        # and mark gene for removal
        self.assertEquals(2, len(self.test_gene1.mrnas))
        self.test_gene1.trim([0, 0])
        self.assertEquals(0, len(self.test_gene1.mrnas))
        self.assertEquals([0, 0], self.test_gene1.indices)

    def test_to_tbl_positive(self):
        gene = Gene(seq_name="seq1", source="maker", indices=[1, 50], strand="+", identifier="foo_gene_1")
        self.assertFalse(gene.annotations)
        gene.add_annotation('foo', 'dog')
        mrna1 = Mock()
        mrna1.to_tbl.return_value = "mrna1_to_tbl...\n"
        mrna2 = Mock()
        mrna2.to_tbl.return_value = "mrna2_to_tbl...\n"
        gene.add_mrna(mrna1)
        gene.add_mrna(mrna2)
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
        gene.add_mrna(mrna1)
        gene.add_mrna(mrna2)
        expected = "50\t1\tgene\n\t\t\tlocus_tag\tfoo_gene_1\nmrna1_to_tbl...\nmrna2_to_tbl...\n"
        self.assertEquals(gene.to_tbl(), expected)

        


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGene))
    return suite

if __name__ == '__main__':
    unittest.main()

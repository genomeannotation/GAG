#!/usr/bin/env python

import unittest
from mock import Mock, PropertyMock
from src.feature_classes import GenePart, CDS, Exon, MRNA, Gene

class TestGene(unittest.TestCase):

    def setUp(self):
        self.test_gene0 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[3734, 7436], strand='+', identifier=1, name="BDOR_007864")
        self.test_gene1 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[3734, 7436], strand='+', identifier=1, name="BDOR_007864")
        self.fake_mrna1 = Mock()
        self.fake_mrna2 = Mock()
        self.test_gene1.add_mrna(self.fake_mrna1)
        self.test_gene1.add_mrna(self.fake_mrna2)


    def test_constructor(self):
        self.assertEqual('Gene', self.test_gene0.__class__.__name__)

    def test_length(self):
        self.assertEqual(3703, self.test_gene0.length())

    def test_add_mrna(self):
        self.assertEquals(0, len(self.test_gene0.mrnas))
        self.test_gene0.add_mrna(self.fake_mrna1)
        self.assertEquals(1, len(self.test_gene0.mrnas))
        self.test_gene0.add_mrna(self.fake_mrna2)
        self.assertEquals(2, len(self.test_gene0.mrnas))

    def test_length_of_shortest_cds_segment(self):
        self.fake_mrna1.length_of_shortest_cds_segment.return_value = 358
        self.fake_mrna2.length_of_shortest_cds_segment.return_value = 241
        self.assertEquals(241, self.test_gene1.length_of_shortest_cds_segment())
        self.fake_mrna1.length_of_shortest_cds_segment.assert_called_with()
        self.fake_mrna2.length_of_shortest_cds_segment.assert_called_with()

    def test_adjust_indices(self):
        self.test_gene1.adjust_indices(16)
        self.fake_mrna1.adjust_indices.assert_called_with(16)
        self.assertEquals(3750, self.test_gene1.indices[0])
        # adjust them back
        self.test_gene1.adjust_indices(-16)
        self.fake_mrna1.adjust_indices.assert_called_with(-16)
        self.assertEquals(3734, self.test_gene1.indices[0])

    def test_collides_range(self):
        self.assertEquals(self.test_gene1.collidesRange(3500, 3800), True)
        self.assertEquals(self.test_gene1.collidesRange(3500, 3600), False)

    def test_to_gff(self):
        self.fake_mrna1.to_gff.return_value = "fake mrna1 to gff here:)\n"
        self.fake_mrna2.to_gff.return_value = "fake mrna2 to gff here:)\n"
        expected = "sctg_0080_0020\tmaker\tgene\t3734\t7436\t.\t+\t."
        expected += "\tID=1;Name=BDOR_007864\n"
        expected += "fake mrna1 to gff here:)\n"
        expected += "fake mrna2 to gff here:)\n"
        self.assertEquals(expected, self.test_gene1.to_gff())

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
        self.fake_mrna1.adjust_indices.assert_called_with(-9)
        self.fake_mrna2.adjust_indices.assert_called_with(-9)
        
    def test_clean_up_indices(self):
        # if indices[0] < 1, set to 1 
        nice_gene = Gene(seq_name='fooseq', source='maker', indices=[-23, 127], strand='-', identifier='foo', name='foo')
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
        junk_gene = Gene(seq_name='barseq', source='maker', indices=[-400, -100], strand='-', identifier='bar', name='bar')
        junk_gene.clean_up_indices()
        self.assertEquals(0, junk_gene.indices[0])
        self.assertEquals(0, junk_gene.indices[1])
        
    def test_remove_invalid_features(self):
        gene = Gene(seq_name="sctg_foo", source='maker', indices=[100, 200], strand='-', identifier='foo_gene', name='gene1')
        nice_mrna = Mock()
        nice_ind = PropertyMock(return_value = [100, 200])
        type(nice_mrna).indices = nice_ind
        junk_mrna = MRNA(identifier='junk', name='junk', indices=[0, 0], parent_id='foo_gene')
        gene.add_mrna(junk_mrna)
        gene.add_mrna(nice_mrna) 
        self.assertEquals(2, len(gene.mrnas))
        # should remove junk_mrna, keeping nice_mrna
        gene.remove_invalid_features()
        self.assertEquals(1, len(gene.mrnas))
        # should call recursively on any mrnas not discarded
        nice_mrna.remove_invalid_features.assert_called_with()

    def test_trim(self):
        gene = Gene(seq_name='sctg_foo', source='maker', indices=[100, 200], strand='-', identifier='foo_gene', name='gene1')
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
        nice_mrna.adjust_indices.assert_called_with(-119)
        bad_mrna.adjust_indices.assert_called_with(-119)
        nice_mrna.adjust_phase.assert_called_with()
        bad_mrna.adjust_phase.assert_called_with()
        nice_mrna.clean_up_indices.assert_called_with()
        bad_mrna.clean_up_indices.assert_called_with()
        nice_mrna.remove_invalid_features.assert_called_with()
        # verify bad_mrna was removed
        self.assertEquals(1, len(gene.mrnas))

        
        

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGene))
    return suite

if __name__ == '__main__':
    unittest.main()

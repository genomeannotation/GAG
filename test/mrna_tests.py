#!/usr/bin/env python

import unittest
from mock import Mock, PropertyMock
from src.gene_part import GenePart, CDS, Exon
from src.mrna import MRNA
from src.gene import Gene

class TestMRNA(unittest.TestCase):

    def setUp(self):
        self.test_mrna0 = MRNA(identifier=2, name="BDOR_007864-RA", indices=[3734, 7436], parent_id=1)
        self.test_mrna1 = MRNA(identifier=2, name="BDOR_007864-RA", indices=[3734, 7436], parent_id=1)
        self.fake_exon = Mock()
        self.fake_cds = Mock()
        self.fake_start_codon = Mock()
        self.test_mrna1.set_exon(self.fake_exon)
        self.test_mrna1.set_cds(self.fake_cds)
        self.test_mrna1.add_other_feature(self.fake_start_codon)

    def test_constructor(self):
        self.assertEqual('MRNA', self.test_mrna0.__class__.__name__)

    def test_length(self):
        self.assertEqual(3703, self.test_mrna0.length())

    def test_set_exon(self):
        self.assertFalse(self.test_mrna0.exon)
        self.test_mrna0.set_exon(self.fake_exon)
        self.assertTrue(self.test_mrna0.exon)

    def test_set_cds(self):
        self.assertFalse(self.test_mrna0.cds)
        self.test_mrna0.set_cds(self.fake_cds)
        self.assertTrue(self.test_mrna0.cds)
       
    def test_add_other_feature(self): 
        self.assertEquals(0, len(self.test_mrna0.other_features))
        self.test_mrna0.add_other_feature(self.fake_start_codon)
        self.assertEquals(1, len(self.test_mrna0.other_features))

    def test_add_start_codon(self):
        self.assertFalse(self.test_mrna0.has_start())
        #self.test_mrna0 = MRNA(identifier=2, name="BDOR_007864-RA", indices=[3734, 7436], parent_id=1)
        self.test_mrna0.add_start_codon(4000)
        self.assertTrue(self.test_mrna0.has_start())
        self.assertEquals([[4000, 4002]], self.test_mrna0.other_features[0].indices)

    def test_add_stop_codon(self):
        self.assertFalse(self.test_mrna0.has_stop())
        self.test_mrna0.add_stop_codon(7002)
        self.assertTrue(self.test_mrna0.has_stop())

    def test_get_cds_indices(self):
        self.fake_cds.indices = [[4000, 4100], [5000, 7000]]
        self.fake_cds.get_phase.return_value = 0
        self.assertEquals([[4000, 4100], [5000, 7000]], self.test_mrna1.get_cds_indices())

    def test_get_cds_when_phase_matters(self):
        self.fake_cds.indices = [[4000, 4100], [5000, 7000]]
        self.fake_cds.get_phase.return_value = 1
        self.assertEquals([[4001, 4100], [5001, 7000]], self.test_mrna1.get_cds_indices())

    def test_adjust_indices(self):
        self.test_mrna1.adjust_indices(32)
        self.assertEqual(7468, self.test_mrna1.indices[1])
        self.fake_cds.adjust_indices.assert_called_with(32)
        self.fake_exon.adjust_indices.assert_called_with(32)
        self.fake_start_codon.adjust_indices.assert_called_with(32)
        # adjust 'em back
        self.test_mrna1.adjust_indices(-32)
        self.assertEqual(3734, self.test_mrna1.indices[0])
        self.fake_cds.adjust_indices.assert_called_with(-32)
        self.fake_exon.adjust_indices.assert_called_with(-32)
        self.fake_start_codon.adjust_indices.assert_called_with(-32)

    def test_length_of_shortest_cds_segment(self):
        self.fake_cds.length_of_shortest_segment.return_value = 241
        self.assertEquals(241, self.test_mrna1.length_of_shortest_cds_segment())
        self.fake_cds.length_of_shortest_segment.assert_called_once_with()

    def test_has_start(self):
        self.fake_start_codon.feature_type = 'start_codon'
        self.assertTrue(self.test_mrna1.has_start())

    def test_has_stop(self):
        self.assertFalse(self.test_mrna1.has_stop())

    def test_str(self):
        expected = "mRNA (ID=2, Name=BDOR_007864-RA) containing Exon, CDS and 1 other features"
        self.assertEquals(expected, str(self.test_mrna1))

    def test_to_gff(self):
        self.fake_exon.to_gff.return_value = "...exon to gff\n"
        self.fake_cds.to_gff.return_value = "...cds to gff\n"
        self.fake_start_codon.to_gff.return_value = "...start codon to gff\n"
        expected = "sctg_0080_0020\tmaker\tmRNA\t"
        expected += "3734\t7436\t.\t+\t.\t"
        expected += "ID=2;Name=BDOR_007864-RA;Parent=1\n"
        expected += "...exon to gff\n...cds to gff\n"
        expected += "...start codon to gff\n"
        actual = self.test_mrna1.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEquals(expected, actual)
        self.fake_exon.to_gff.assert_called_with("sctg_0080_0020", "maker", '+')
        self.fake_cds.to_gff.assert_called_with("sctg_0080_0020", "maker", '+')
        self.fake_start_codon.to_gff.assert_called_with("sctg_0080_0020", "maker", '+')

    def test_remove_first_cds_segment_if_shorter_than(self):
        bad_indices = [[5,7], [10,20]]
        self.fake_cds.indices = bad_indices
        self.assertEquals(2, len(self.test_mrna1.cds.indices))
        self.test_mrna1.remove_first_cds_segment_if_shorter_than(4)
        self.assertEquals(1, len(self.test_mrna1.cds.indices))

    def test_remove_first_cds_segment_if_shorter_than_leaves_well_enough_alone(self):
        good_indices = [[5, 8], [10, 20]]
        self.fake_cds.indices = good_indices
        self.assertEquals(2, len(self.test_mrna1.cds.indices))
        self.test_mrna1.remove_first_cds_segment_if_shorter_than(4)
        self.assertEquals(2, len(self.test_mrna1.cds.indices))

    def test_trim_end(self):
        self.test_mrna1.trim_end(7400)
        self.assertEquals(3734, self.test_mrna1.indices[0])
        self.assertEquals(7400, self.test_mrna1.indices[1])
        # verify recursive call on child features
        self.fake_exon.trim_end.assert_called_with(7400)
        self.fake_cds.trim_end.assert_called_with(7400)
        self.fake_start_codon.trim_end.assert_called_with(7400)
        # if trim before mrna begin, mark for removal
        #    by changing indices to [0, 0]
        self.test_mrna1.trim_end(50)
        self.assertEquals(0, self.test_mrna1.indices[0])
        self.assertEquals(0, self.test_mrna1.indices[1])

    def test_adjust_phase(self):
        # just need to verify call to child cds
        self.test_mrna1.adjust_phase()
        self.fake_cds.adjust_phase.assert_called_with()

    def test_clean_up_indices(self):
        nice_mrna = MRNA(identifier='foo', name='foo', indices=[-10, 200], parent_id='foo')
        junk_mrna = MRNA(identifier='bar', name='bar', indices=[-300, -200], parent_id='bar')
        nice_mrna.clean_up_indices()
        self.assertEquals(1, nice_mrna.indices[0])
        self.assertEquals(200, nice_mrna.indices[1])
        junk_mrna.clean_up_indices()
        self.assertEquals(0, junk_mrna.indices[0])
        self.assertEquals(0, junk_mrna.indices[1])
        # test recursive call
        codon = Mock()
        exon = Mock()
        nice_mrna.add_other_feature(codon)
        nice_mrna.set_exon(exon)
        nice_mrna.clean_up_indices()
        codon.clean_up_indices.assert_called_with()
        exon.clean_up_indices.assert_called_with()

    def test_create_start_and_stop_if_necessary(self):
        fasta = Mock()
        seq_name = 'seq_foo'
        strand = '+'
        self.test_mrna0.create_start_and_stop_if_necessary(fasta, seq_name, strand)
        pass
        #self.fake_cds = Mock()
        #self.test_mrna0 = MRNA(identifier=2, name="BDOR_007864-RA", indices=[3734, 7436], parent_id=1)
            #mrna.create_start_and_stop_if_necessary(fasta, self.strand)

    def test_remove_invalid_features(self):
        mrna = MRNA(identifier='foo', name='foo', indices=[20, 50], parent_id='foo') 
        # set up mocks...
        nice_exon = Mock()
        nice1 = PropertyMock(return_value = [10, 50])
        type(nice_exon).indices = nice1

        empty_cds = Mock() 
        empty_stop_codon = Mock()
        empty = PropertyMock(return_value=[])
        type(empty_cds).indices = empty
        type(empty_stop_codon).indices = empty

        invalid_start_codon = Mock()
        invalid_start_codon.valid_codon.return_value = False
        invalid1 = PropertyMock(return_value = [1, 2])
        starttype = PropertyMock(return_value = 'start_codon')
        type(invalid_start_codon).indices = invalid1
        type(invalid_start_codon).feature_type = starttype

        nice_utr = Mock()
        nice2 = PropertyMock(return_value = [1, 9])
        type(nice_utr).indices = nice2

        # add mock features
        mrna.set_exon(nice_exon)
        mrna.set_cds(empty_cds)
        mrna.add_other_feature(empty_stop_codon)
        mrna.add_other_feature(nice_utr)
        mrna.add_other_feature(invalid_start_codon)

        # (finally) run test
        self.assertEquals(3, len(mrna.other_features))
        self.assertTrue(mrna.cds)
        self.assertTrue(mrna.exon)
        mrna.remove_invalid_features()
        self.assertFalse(mrna.cds)
        self.assertTrue(mrna.exon)
        self.assertEquals(1, len(mrna.other_features))

    def test_is_maker_mrna(self):
        maker_mrna = MRNA(identifier=1, name="maker-scaffold00080-est_gff_Cufflinks-gene-2.9-mRNA-1", indices=[3734, 7436], parent_id=1)
        self.assertFalse(self.test_mrna0.is_maker_mrna())
        self.assertTrue(maker_mrna.is_maker_mrna())
       



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMRNA))
    return suite

if __name__ == '__main__':
    unittest.main()

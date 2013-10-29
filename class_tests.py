#!/usr/bin/env python

import unittest
from mock import Mock
from feature_classes import CDS, Exon, OtherFeature, MRNA, Gene

class TestFeatureClasses(unittest.TestCase):

    def test_CDS(self):
        # test constructor
        test_indices1 = [[3734, 4034], [4092, 4332], [4399, 5185], [5249, 6565], [6630, 7436]]
        test_names1 = ["BDOR_007864-RA:cds:0", "BDOR_007864-RA:cds:1", "BDOR_007864-RA:cds:2", "BDOR_007864-RA:cds:3", "BDOR_007864-RA:cds:4"]
        test_ids1 = [8, 9, 10, 11, 12]
        test_frames1 = [0, 2, 1, 0, 0]
        test_parent_id1 = 2
        test_cds1 = CDS(ids=test_ids1, names=test_names1, indices=test_indices1, frames=test_frames1, parent_id=test_parent_id1)
        self.assertEqual('CDS', test_cds1.__class__.__name__)

        # test .length_of_shortest_segment
        self.assertEqual(241, test_cds1.length_of_shortest_segment())

        # test .length
        self.assertEqual(3453, test_cds1.length())

        # test .adjust_indices
        test_cds1.adjust_indices(-5)
        self.assertEqual(3729, test_cds1.indices[0][0])
        # (adjust them back so future test don't get confused :)
        test_cds1.adjust_indices(5)
        self.assertEqual(5185, test_cds1.indices[2][1])

        # test to_gff
        expected1 = "sctg_0080_0020\tmaker\tCDS\t3734\t4034\t.\t+\t0\tID=8;Name=BDOR_007864-RA:cds:0;Parent=2\n"
        expected2 = "sctg_0080_0020\tmaker\tCDS\t4092\t4332\t.\t+\t2\tID=9;Name=BDOR_007864-RA:cds:1;Parent=2\n"
        expected3 = "sctg_0080_0020\tmaker\tCDS\t4399\t5185\t.\t+\t1\tID=10;Name=BDOR_007864-RA:cds:2;Parent=2\n"
        expected4 = "sctg_0080_0020\tmaker\tCDS\t5249\t6565\t.\t+\t0\tID=11;Name=BDOR_007864-RA:cds:3;Parent=2\n"
        expected5 = "sctg_0080_0020\tmaker\tCDS\t6630\t7436\t.\t+\t0\tID=12;Name=BDOR_007864-RA:cds:4;Parent=2\n"
        expected = expected1 + expected2 + expected3 + expected4 + expected5
        actual = test_cds1.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEqual(expected, actual)
        # what if id, parent_id are strings? does it matter?
        test_cds2 = CDS(ids=['foo1', 'foo2', 'foo3', 'foo4', 'foo5'], names=test_names1, indices=test_indices1, frames=test_frames1, parent_id='bar7')
        expected6 = "sctg_0080_0020\tmaker\tCDS\t3734\t4034\t.\t+\t0\tID=foo1;Name=BDOR_007864-RA:cds:0;Parent=bar7\n"
        expected7 = "sctg_0080_0020\tmaker\tCDS\t4092\t4332\t.\t+\t2\tID=foo2;Name=BDOR_007864-RA:cds:1;Parent=bar7\n"
        expected8 = "sctg_0080_0020\tmaker\tCDS\t4399\t5185\t.\t+\t1\tID=foo3;Name=BDOR_007864-RA:cds:2;Parent=bar7\n"
        expected9 = "sctg_0080_0020\tmaker\tCDS\t5249\t6565\t.\t+\t0\tID=foo4;Name=BDOR_007864-RA:cds:3;Parent=bar7\n"
        expected10 = "sctg_0080_0020\tmaker\tCDS\t6630\t7436\t.\t+\t0\tID=foo5;Name=BDOR_007864-RA:cds:4;Parent=bar7\n"
        expected = expected6 + expected7 + expected8 + expected9 + expected10
        actual = test_cds2.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEqual(expected, actual)

        # test to_tbl
        # TODO when add annotations stuff

    def test_Exon(self):
        # test constructor
        test_ids1 = [3, 4, 5, 6, 7]
        test_names1 = ["BDOR_007864-RA:exon:0", "BDOR_007864-RA:exon:1", "BDOR_007864-RA:exon:2", "BDOR_007864-RA:exon:3", "BDOR_007864-RA:exon:4"]
        test_indices1 = [[3734, 4034], [4092, 4332], [4399, 5185], [5249, 6565], [6630, 7436]]
        test_scores1 = [0.9, 0.9, 0.9, 0.9, 0.9]
        test_parent_id1 = 2
        test_exon1 = Exon(ids=test_ids1, names=test_names1, indices=test_indices1, scores=test_scores1, parent_id=test_parent_id1)
        self.assertEqual('Exon', test_exon1.__class__.__name__)
        
        # test .length
        self.assertEqual(3453, test_exon1.length())

        # test .adjust_indices
        test_exon1.adjust_indices(7)
        self.assertEqual(4339, test_exon1.indices[1][1])
        # adjust them back...
        test_exon1.adjust_indices(-7)
        self.assertEqual(4332, test_exon1.indices[1][1])

        # test .to_gff
        expected1 = "sctg_0080_0020\tmaker\texon\t3734\t4034\t0.9\t+\t.\tID=3;Name=BDOR_007864-RA:exon:0;Parent=2\n"
        expected2 = "sctg_0080_0020\tmaker\texon\t4092\t4332\t0.9\t+\t.\tID=4;Name=BDOR_007864-RA:exon:1;Parent=2\n"
        expected3 = "sctg_0080_0020\tmaker\texon\t4399\t5185\t0.9\t+\t.\tID=5;Name=BDOR_007864-RA:exon:2;Parent=2\n"
        expected4 = "sctg_0080_0020\tmaker\texon\t5249\t6565\t0.9\t+\t.\tID=6;Name=BDOR_007864-RA:exon:3;Parent=2\n"
        expected5 = "sctg_0080_0020\tmaker\texon\t6630\t7436\t0.9\t+\t.\tID=7;Name=BDOR_007864-RA:exon:4;Parent=2\n"
        expected = expected1 + expected2 + expected3 + expected4 + expected5
        actual = test_exon1.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEqual(expected, actual)


    def test_OtherFeature(self):
        # test constructor
        test_start1 = OtherFeature(feature_type="start_codon", indices=[3734, 3736], id=13, name="BDOR_007864-RA:start1", parent_id=2)
        self.assertEqual('OtherFeature', test_start1.__class__.__name__)

        # test .length
        self.assertEqual(3, test_start1.length())

        # test .adjust_indices
        test_start1.adjust_indices(-20)
        self.assertEqual(3714, test_start1.indices[0])
        # adjust them back for consistency...
        test_start1.adjust_indices(20)
        self.assertEqual(3736, test_start1.indices[1])

        # test .to_gff
        expected = "sctg_0080_0020\tmaker\tstart_codon\t"
        expected += "3734\t3736\t.\t+\t."
        expected += "ID=13;Name=BDOR_007864-RA:start1;Parent=2\n"
        actual = test_start1.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEqual(expected, actual)

    def test_MRNA(self):
        # test constructor
        test_mrna1 = MRNA(id=2, name="BDOR_007864-RA", indices=[3734, 7436], parent_id=1)
        self.assertEqual('MRNA', test_mrna1.__class__.__name__)

        # test .length
        self.assertEqual(3703, test_mrna1.length())

        # test .adjust_indices
        test_mrna1.adjust_indices(32)
        self.assertEqual(7468, test_mrna1.indices[1])
        # adjust 'em back
        test_mrna1.adjust_indices(-32)
        self.assertEqual(3734, test_mrna1.indices[0])

        # TODO should call adjust indices on exon, cds and other features!

        # test .set_exon
        self.assertFalse(test_mrna1.exon)
        fake_exon = Mock()
        test_mrna1.set_exon(fake_exon)
        self.assertTrue(test_mrna1.exon)

        # test .set_cds
        self.assertFalse(test_mrna1.cds)
        fake_cds = Mock()
        test_mrna1.set_cds(fake_cds)
        self.assertTrue(test_mrna1.cds)

        # test .add_other_feature
        self.assertEquals(0, len(test_mrna1.other_features))
        fake_start_codon = Mock()
        test_mrna1.add_other_feature(fake_start_codon)
        self.assertEquals(1, len(test_mrna1.other_features))
        fake_stop_codon = Mock()
        test_mrna1.add_other_feature(fake_stop_codon)
        self.assertEquals(2, len(test_mrna1.other_features))

        # test .length_of_shortest_cds_segment
        fake_cds.length_of_shortest_segment.return_value = 241
        self.assertEquals(241, test_mrna1.length_of_shortest_cds_segment())
        fake_cds.length_of_shortest_segment.assert_called()

        # test .to_gff (this is where the money is)
        fake_exon.to_gff.return_value = "...exon to gff\n"
        fake_cds.to_gff.return_value = "...cds to gff\n"
        fake_start_codon.to_gff.return_value = "...start codon to gff\n"
        fake_stop_codon.to_gff.return_value = "...stop codon to gff\n"
        expected = "sctg_0080_0020\tmaker\tmRNA\t"
        expected += "3734\t7436\t.\t+\t.\t"
        expected += "ID=2;Name=BDOR_007864-RA;Parent=1\n"
        expected += "...exon to gff\n...cds to gff\n"
        expected += "...start codon to gff\n...stop codon to gff\n"
        actual = test_mrna1.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEquals(expected, actual)
        fake_exon.to_gff.assert_called_with("sctg_0080_0020", "maker", '+')
        fake_cds.to_gff.assert_called_with("sctg_0080_0020", "maker", '+')
        fake_start_codon.to_gff.assert_called_with("sctg_0080_0020", "maker", '+')
        fake_stop_codon.to_gff.assert_called_with("sctg_0080_0020", "maker", '+')


    def test_Gene(self):
        # test constructor
        test_gene1 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[3734, 7436], strand='+', id=1, name="BDOR_007864")
        self.assertEqual('Gene', test_gene1.__class__.__name__)

        # test .length
        self.assertEqual(3703, test_gene1.length())

        # test .add_mrna
        self.assertEquals(0, len(test_gene1.mrnas))
        fake_mrna1 = Mock()
        test_gene1.add_mrna(fake_mrna1)
        self.assertEquals(1, len(test_gene1.mrnas))
        fake_mrna2 = Mock()
        test_gene1.add_mrna(fake_mrna2)
        self.assertEquals(2, len(test_gene1.mrnas))

        # test .length_of_shortest_cds_segment
        fake_mrna1.length_of_shortest_cds_segment.return_value = 358
        fake_mrna2.length_of_shortest_cds_segment.return_value = 241
        self.assertEquals(241, test_gene1.length_of_shortest_cds_segment())
        fake_mrna1.length_of_shortest_cds_segment.assert_called()
        fake_mrna2.length_of_shortest_cds_segment.assert_called()

        # TODO test adjust indices (must call recursively!)
        # TODO test has_stop, has_start (must add to MRNA first)
        # TODO to_gff
        
        
        
        
        
        

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFeatureClasses))
    return suite

if __name__ == '__main__':
    unittest.main()

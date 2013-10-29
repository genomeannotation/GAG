#!/usr/bin/env python

import unittest
from mock import Mock
from feature_classes import GenePart, CDS, Exon, MRNA, Gene

class TestFeatureClasses(unittest.TestCase):

    def test_GenePart(self):
        # test constructor
        gp1 = GenePart()
        self.assertEqual(0, len(gp1.indices))
        self.assertFalse(gp1.score)
        self.assertFalse(gp1.parent_id)
        gp2 = GenePart(feature_type='CDS', indices=[[1, 44], [65, 103]])
        self.assertEqual(2, len(gp2.indices))
        self.assertEqual('CDS', gp2.feature_type)

        # test .length
        self.assertEqual(83, gp2.length())
        # what if no indices at all?
        self.assertFalse(gp1.length())

        # test .adjust_indices
        gp2.adjust_indices(10)
        self.assertEqual(54, gp2.indices[0][1])
        # now put it back
        gp2.adjust_indices(-10)
        self.assertEqual(65, gp2.indices[1][0])

        # test .length_of_shortest_segment
        self.assertEqual(39, gp2.length_of_shortest_segment())
        # what if no indices?
        self.assertFalse(gp1.length_of_shortest_segment())

        # test .generate_attribute_entry
        gp2.id = ['foo1', 'foo2']
        gp2.parent_id = 'mama'
        expected = "ID=foo2;Parent=mama\n"
        self.assertEqual(expected, gp2.generate_attribute_entry(1))
        gp2.name = ['cds:0', 'cds:1']
        expected = "ID=foo1;Name=cds:0;Parent=mama\n"
        self.assertEqual(expected, gp2.generate_attribute_entry(0))
        # what if index out of range?
        self.assertFalse(gp2.generate_attribute_entry(2))
        # what if no id or parent_id?
        self.assertFalse(gp1.generate_attribute_entry(0))
        gp1.parent_id = 'dad'
        self.assertFalse(gp1.generate_attribute_entry(0))

        # test .to_gff
        seq_name = "sctg_0001_0001"
        source = 'maker'
        strand = '+'
        expected = "sctg_0001_0001\tmaker\tCDS\t1\t44\t.\t+\t.\t"
        expected += "ID=foo1;Name=cds:0;Parent=mama\n"
        expected += "sctg_0001_0001\tmaker\tCDS\t65\t103\t.\t+\t.\t"
        expected += "ID=foo2;Name=cds:1;Parent=mama\n"
        actual = gp2.to_gff(seq_name=seq_name, source=source, strand=strand)
        self.assertEqual(expected, actual)
        # what if no indices, etc.?
        self.assertFalse(gp1.to_gff(seq_name="foo", source="bar", strand=":)"))
        

    def test_CDS(self):
        # test constructor
        test_indices1 = [[3734, 4034], [4092, 4332], [4399, 5185], [5249, 6565], [6630, 7436]]
        test_name1 = ["BDOR_007864-RA:cds:0", "BDOR_007864-RA:cds:1", "BDOR_007864-RA:cds:2", "BDOR_007864-RA:cds:3", "BDOR_007864-RA:cds:4"]
        test_id1 = [8, 9, 10, 11, 12]
        test_phase1 = [0, 2, 1, 0, 0]
        test_parent_id1 = 2
        test_cds1 = CDS(id=test_id1, name=test_name1, indices=test_indices1, score=None, phase=test_phase1, parent_id=test_parent_id1)
        self.assertEqual('CDS', test_cds1.__class__.__name__)
        # should also be able to construct w/o all the params...
        empty_cds = CDS()
        self.assertEqual('CDS', empty_cds.feature_type) 

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
        test_cds2 = CDS(id=['foo1', 'foo2', 'foo3', 'foo4', 'foo5'], name=test_name1, indices=test_indices1, score=None, phase=test_phase1, parent_id='bar7')
        expected1 = "sctg_0080_0020\tmaker\tCDS\t3734\t4034\t.\t+\t0\tID=foo1;Name=BDOR_007864-RA:cds:0;Parent=bar7\n"
        expected2 = "sctg_0080_0020\tmaker\tCDS\t4092\t4332\t.\t+\t2\tID=foo2;Name=BDOR_007864-RA:cds:1;Parent=bar7\n"
        expected3 = "sctg_0080_0020\tmaker\tCDS\t4399\t5185\t.\t+\t1\tID=foo3;Name=BDOR_007864-RA:cds:2;Parent=bar7\n"
        expected4 = "sctg_0080_0020\tmaker\tCDS\t5249\t6565\t.\t+\t0\tID=foo4;Name=BDOR_007864-RA:cds:3;Parent=bar7\n"
        expected5 = "sctg_0080_0020\tmaker\tCDS\t6630\t7436\t.\t+\t0\tID=foo5;Name=BDOR_007864-RA:cds:4;Parent=bar7\n"
        expected = expected1 + expected2 + expected3 + expected4 + expected5
        actual = test_cds2.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEqual(expected, actual)
        # what if name=None?
        test_cds2.name=[]
        expected1 = "sctg_0080_0020\tmaker\tCDS\t3734\t4034\t.\t+\t0\tID=foo1;Parent=bar7\n"
        expected2 = "sctg_0080_0020\tmaker\tCDS\t4092\t4332\t.\t+\t2\tID=foo2;Parent=bar7\n"
        expected3 = "sctg_0080_0020\tmaker\tCDS\t4399\t5185\t.\t+\t1\tID=foo3;Parent=bar7\n"
        expected4 = "sctg_0080_0020\tmaker\tCDS\t5249\t6565\t.\t+\t0\tID=foo4;Parent=bar7\n"
        expected5 = "sctg_0080_0020\tmaker\tCDS\t6630\t7436\t.\t+\t0\tID=foo5;Parent=bar7\n"
        expected = expected1 + expected2 + expected3 + expected4 + expected5
        actual = test_cds2.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEqual(expected, actual)

        # test to_tbl
        # TODO when add annotations stuff

    def test_Exon(self):
        # test constructor
        test_id1 = [3, 4, 5, 6, 7]
        test_name1 = ["BDOR_007864-RA:exon:0", "BDOR_007864-RA:exon:1", "BDOR_007864-RA:exon:2", "BDOR_007864-RA:exon:3", "BDOR_007864-RA:exon:4"]
        test_indices1 = [[3734, 4034], [4092, 4332], [4399, 5185], [5249, 6565], [6630, 7436]]
        test_score1 = [0.9, 0.9, 0.9, 0.9, 0.9]
        test_parent_id1 = 2
        test_exon1 = Exon(id=test_id1, name=test_name1, indices=test_indices1, score=test_score1, parent_id=test_parent_id1)
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
        # TODO Gene.adjust_indices should be robust against attempts to 
        #  adjust into negative territory; it should modify downstream calls
        #  in order to compensate. hence no error checking for this kind of thing
        #  in CDS.adjust-indices; must trust call from Gene.
        
        
        
        
        
        

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFeatureClasses))
    return suite

if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python

import unittest
from mock import Mock, PropertyMock
from src.gene_part import GenePart, CDS, Exon
from src.mrna import MRNA
from src.gene import Gene

class TestGenePart(unittest.TestCase):

    def setUp(self):
        self.gp1 = GenePart()
        self.gp2 = GenePart(feature_type='CDS', indices=[1, 44])
        self.gp2.add_indices([65, 103])
        self.gp2.name = ['cds:0', 'cds:1']
        self.gp2.identifier = ['foo1', 'foo2']
        self.gp2.parent_id = 'mama'
        self.gp3 = GenePart(feature_type='exon')

    def test_constructor(self):
        self.gp1 = GenePart()
        self.assertEquals(0, len(self.gp1.indices))
        self.assertFalse(self.gp1.score)
        self.assertFalse(self.gp1.parent_id)
        self.assertEquals('CDS', self.gp2.feature_type)

    def test_add_indices(self):
        self.assertEquals(2, len(self.gp2.indices))
        self.gp2.add_indices([150, 197])
        self.assertEquals(3, len(self.gp2.indices))

        self.assertEquals(0, len(self.gp3.indices))
        self.gp3.add_indices([77, 144])
        self.assertEquals(1, len(self.gp3.indices))
        # error check
        self.assertRaises(ValueError, self.gp3.add_indices, 7)
        self.assertRaises(ValueError, self.gp3.add_indices, 'foo')

    def test_add_name(self):
        self.assertEquals(0, len(self.gp3.name))
        self.gp3.add_name('BDOR_007864-RA:cds:0')
        self.assertEquals(1, len(self.gp3.name))

    def test_add_identifier(self):
        self.assertEquals(0, len(self.gp3.identifier))
        self.gp3.add_identifier('7')
        self.assertEquals(1, len(self.gp3.identifier))

    def test_length(self):
        self.assertEquals(83, self.gp2.length())
        # what if no indices at all?
        self.assertFalse(self.gp1.length())

    def test_valid_codon(self):
        bad_stop = GenePart(feature_type='stop_codon', identifier='foo', name='foo', indices=[2, 3], parent_id='foofoo')
        self.assertFalse(bad_stop.valid_codon())
        good_stop = GenePart(feature_type='stop_codon', identifier='foo', name='foo', indices=[2, 4], parent_id='foofoo')
        self.assertTrue(good_stop.valid_codon())

    def test_adjust_indices(self):
        self.gp2.adjust_indices(10)
        self.assertEquals(54, self.gp2.indices[0][1])
        # now put it back
        self.gp2.adjust_indices(-10)
        self.assertEquals(65, self.gp2.indices[1][0])

    def test_adjust_indices_after_start_index(self):
        self.gp2.adjust_indices(10, 50)
        self.assertEqual(44, self.gp2.indices[0][1])
        self.assertEqual(75, self.gp2.indices[1][0])

    def test_clean_up_indices(self):
        # test on a GenePart...
        nice_cds = CDS(identifier='foo', name='foo', indices=[-5, 28], phase=1, parent_id='bar')
        junk_cds = CDS(identifier='bar', name='bar', indices=[-20, -2], phase=0, parent_id='bar')
        nice_cds.clean_up_indices()
        self.assertEquals(1, nice_cds.indices[0][0])
        self.assertEquals(28, nice_cds.indices[0][1]) 
        junk_cds.clean_up_indices()
        self.assertEquals(0, junk_cds.indices[0][0])
        self.assertEquals(0, junk_cds.indices[0][1])

    def test_length_of_shortest_segment(self):
        self.assertEquals(39, self.gp2.length_of_shortest_segment())
        # what if no indices?
        self.assertFalse(self.gp1.length_of_shortest_segment())

    def test_trim_end(self):
        self.gp2.trim_end(100)
        self.assertEquals(65, self.gp2.indices[1][0]) #unchanged
        self.assertEquals(100, self.gp2.indices[1][1])
        self.gp2.trim_end(30)
        # should mark second segment for removal
        # by setting indices = [0, 0]
        self.assertEquals([0,0], self.gp2.indices[1])
        self.assertEquals([1, 30], self.gp2.indices[0])

    def test_remove_trimmed_segments(self):
        testgp = Exon(identifier='foo', name='foo', indices=[0, 0], parent_id='foo')
        testgp.add_identifier('foo2')
        testgp.add_name('foo2')
        testgp.add_indices([22, 100])
        self.assertEquals(2, len(testgp.identifier))
        self.assertEquals(2, len(testgp.indices))
        testgp.remove_trimmed_segments()
        self.assertEquals(1, len(testgp.name))
        self.assertEquals(1, len(testgp.indices))
        # verify that it will remove all segments if appropriate
        testgp.indices[0] = [0, 0]
        testgp.remove_trimmed_segments()
        self.assertEquals(0, len(testgp.indices))
        self.assertEquals(0, len(testgp.identifier))

    def test_invalidate_region(self):
        expected = [[0, 0], [71, 89]]
        self.gp2.invalidate_region(1, 44)
        self.gp2.invalidate_region(50, 70)
        self.gp2.invalidate_region(90, 110)
        self.assertEquals(expected, self.gp2.indices)

    def test_invalidate_region_chops_off_beginning(self):
        expected = [[5, 44], [65, 103]]
        self.gp2.invalidate_region(1, 4)
        self.assertEqual(expected, self.gp2.indices)

    def test_generate_attribute_entry(self):
        # test .generate_attribute_entry
        expected = "ID=foo1;Name=cds:0;Parent=mama\n"
        self.assertEquals(expected, self.gp2.generate_attribute_entry(0))
        expected = "ID=foo2;Name=cds:1;Parent=mama\n"
        self.assertEquals(expected, self.gp2.generate_attribute_entry(1))
        # what if index out of range?
        self.assertFalse(self.gp2.generate_attribute_entry(2))
        # what if no identifier or parent_id?
        self.assertFalse(self.gp1.generate_attribute_entry(0))
        self.gp1.parent_id = 'dad'
        self.assertFalse(self.gp1.generate_attribute_entry(0))

    def test_str(self):
        expected = "CDS (first ID=foo1, first name=cds:0)"
        self.assertEquals(expected, str(self.gp2))

    def test_to_gff(self):
        # test .to_gff
        seq_name = "sctg_0001_0001"
        source = 'maker'
        strand = '+'
        expected = "sctg_0001_0001\tmaker\tCDS\t1\t44\t.\t+\t.\t"
        expected += "ID=foo1;Name=cds:0;Parent=mama\n"
        expected += "sctg_0001_0001\tmaker\tCDS\t65\t103\t.\t+\t.\t"
        expected += "ID=foo2;Name=cds:1;Parent=mama\n"
        actual = self.gp2.to_gff(seq_name=seq_name, source=source, strand=strand)
        self.assertEquals(expected, actual)
        # what if no indices, etc.?
        self.assertFalse(self.gp1.to_gff(seq_name="foo", source="bar", strand=":)"))
        
class TestCDS(unittest.TestCase):

    def setUp(self):
        self.test_indices1 = [3734, 4034]
        self.extra_indices = [[4092, 4332], [4399, 5185], [5249, 6565], [6630, 7436]]
        self.test_name1 = "BDOR_007864-RA:cds:0"
        self.extra_names = ["BDOR_007864-RA:cds:1", "BDOR_007864-RA:cds:2", "BDOR_007864-RA:cds:3", "BDOR_007864-RA:cds:4"]
        test_identifier1 = 8
        self.extra_identifiers = [9, 10, 11, 12]
        self.test_phase1 = 0
        self.extra_phases = [2, 1, 0, 0]
        test_parent_id1 = 2
        self.test_cds0 = CDS(identifier=test_identifier1, name=self.test_name1, indices=self.test_indices1, score=None, phase=self.test_phase1, parent_id=test_parent_id1)
        self.test_cds1 = CDS(identifier=test_identifier1, name=self.test_name1, indices=self.test_indices1, score=None, phase=self.test_phase1, parent_id=test_parent_id1)
        for ind_pair in self.extra_indices:
            self.test_cds1.add_indices(ind_pair)
        for ident in self.extra_identifiers:
            self.test_cds1.add_identifier(ident)
        for name in self.extra_names:
            self.test_cds1.add_name(name)
        for phase in self.extra_phases:
            self.test_cds1.add_phase(phase)

    def test_get_start_indices_pos_strand(self):
        expected = [3734, 3736]
        self.assertEquals(expected, self.test_cds1.get_start_indices('+'))

    def test_get_start_indices_neg_strand(self):
        expected = [4032, 4034]
        self.assertEquals(expected, self.test_cds1.get_start_indices('-'))

    def test_get_stop_indices_pos_strand(self):
        expected = [7434, 7436]
        self.assertEquals(expected, self.test_cds1.get_stop_indices('+'))

    def test_get_stop_indices_neg_strand(self):
        expected = [6630, 6632]
        self.assertEquals(expected, self.test_cds1.get_stop_indices('-'))

    def test_extract_sequence_pos_strand(self):
        fasta = Mock()
        fasta.get_subseq.return_value = 'GATTACA'
        strand = '+'
        seq_name = 'seq1'
        seq = self.test_cds1.extract_sequence(fasta, seq_name, strand)
        expected = 'GATTACAGATTACAGATTACAGATTACAGATTACA'
        self.assertEquals(expected, seq)

    def test_extract_sequence_neg_strand(self):
        fasta = Mock()
        fasta.get_subseq.return_value = 'GATTACA'
        strand = '-'
        seq_name = 'seq1'
        seq = self.test_cds1.extract_sequence(fasta, seq_name, strand)
        expected = 'TGTAATCTGTAATCTGTAATCTGTAATCTGTAATC'
        self.assertEquals(expected, seq)
        calls = fasta.mock_calls
        # build a list of arguments to all calls to fasta
        first_call_args = calls[0][1]
        self.assertTrue('seq1' in first_call_args)
        second_call_args = calls[1][1]
        self.assertTrue([[4092, 4332]] in second_call_args) 

    def test_cds_constructor(self):
        self.assertEquals('CDS', self.test_cds0.__class__.__name__)
        # should also be able to construct w/o all the params...
        empty_cds = CDS()
        self.assertEquals('CDS', empty_cds.feature_type) 

    def test_add_indices(self):
        for ind_pair in self.extra_indices:
            self.test_cds0.add_indices(ind_pair)
        self.assertEquals([4399, 5185], self.test_cds0.indices[2])

    def test_add_name(self):
        # test .add_name
        for name in self.extra_names:
            self.test_cds0.add_name(name)
        self.assertEquals("BDOR_007864-RA:cds:4", self.test_cds0.name[4])

    def test_add_identifier(self):
        for ident in self.extra_identifiers:
            self.test_cds0.add_identifier(ident)
        self.assertEquals(5, len(self.test_cds0.identifier))

    def test_add_phase(self):
        for phase in self.extra_phases:
            self.test_cds0.add_phase(phase)
        self.assertEquals(5, len(self.test_cds0.phase))
        self.assertEquals(1, self.test_cds0.phase[2])

    def test_adjust_phase(self):
        # formula for phase is (old_phase + negative_index -1) mod 3
        cds1 = CDS(identifier='foo', name='foo_cds', indices=[-4, 13], phase=0, parent_id='bar')
        cds1.adjust_phase()
        self.assertEquals(1, cds1.phase[0])
        cds1.indices[0][0] = -4
        cds1.adjust_phase()
        self.assertEquals(2, cds1.phase[0])
        # shouldn't affect cds with indices[i][0] > 0
        cds1.indices[0] = [1, 10]
        cds1.adjust_phase()
        self.assertEquals([1, 10], cds1.indices[0])

    def test_length_of_shortest_segment(self):
        self.assertEquals(241, self.test_cds1.length_of_shortest_segment())

    def test_length(self):
        self.assertEquals(3453, self.test_cds1.length())

    def test_adjust_indices(self):
        self.test_cds1.adjust_indices(-5)
        self.assertEquals(3729, self.test_cds1.indices[0][0])
        # (adjust them back so future test don't get confused :)
        self.test_cds1.adjust_indices(5)
        self.assertEquals(5185, self.test_cds1.indices[2][1])

    def test_invalidate_region(self):
        expected = [3737, 4034]
        expectedPhase = 2
        self.test_cds1.invalidate_region(3700, 3736)
        self.assertEquals(expected, self.test_cds1.indices[0])
        self.assertEquals(expectedPhase, self.test_cds1.phase[0])

    def test_to_gff(self):
        expected1 = "sctg_0080_0020\tmaker\tCDS\t3734\t4034\t.\t+\t0\tID=8;Name=BDOR_007864-RA:cds:0;Parent=2\n"
        expected2 = "sctg_0080_0020\tmaker\tCDS\t4092\t4332\t.\t+\t2\tID=9;Name=BDOR_007864-RA:cds:1;Parent=2\n"
        expected3 = "sctg_0080_0020\tmaker\tCDS\t4399\t5185\t.\t+\t1\tID=10;Name=BDOR_007864-RA:cds:2;Parent=2\n"
        expected4 = "sctg_0080_0020\tmaker\tCDS\t5249\t6565\t.\t+\t0\tID=11;Name=BDOR_007864-RA:cds:3;Parent=2\n"
        expected5 = "sctg_0080_0020\tmaker\tCDS\t6630\t7436\t.\t+\t0\tID=12;Name=BDOR_007864-RA:cds:4;Parent=2\n"
        expected = expected1 + expected2 + expected3 + expected4 + expected5
        actual = self.test_cds1.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEquals(expected, actual)
        # what if identifier, parent_id are strings? does it matter?
        test_cds2 = CDS(identifier='foo1', name=self.test_name1, indices=self.test_indices1, score=None, phase=self.test_phase1, parent_id='bar7')
        extra_identifiers2 = ['foo2', 'foo3', 'foo4', 'foo5']
        for ind_pair in self.extra_indices:
            test_cds2.add_indices(ind_pair)
        for name in self.extra_names:
            test_cds2.add_name(name)
        for ident in extra_identifiers2:
            test_cds2.add_identifier(ident)
        for phase in self.extra_phases:
            test_cds2.add_phase(phase)
        expected1 = "sctg_0080_0020\tmaker\tCDS\t3734\t4034\t.\t+\t0\tID=foo1;Name=BDOR_007864-RA:cds:0;Parent=bar7\n"
        expected2 = "sctg_0080_0020\tmaker\tCDS\t4092\t4332\t.\t+\t2\tID=foo2;Name=BDOR_007864-RA:cds:1;Parent=bar7\n"
        expected3 = "sctg_0080_0020\tmaker\tCDS\t4399\t5185\t.\t+\t1\tID=foo3;Name=BDOR_007864-RA:cds:2;Parent=bar7\n"
        expected4 = "sctg_0080_0020\tmaker\tCDS\t5249\t6565\t.\t+\t0\tID=foo4;Name=BDOR_007864-RA:cds:3;Parent=bar7\n"
        expected5 = "sctg_0080_0020\tmaker\tCDS\t6630\t7436\t.\t+\t0\tID=foo5;Name=BDOR_007864-RA:cds:4;Parent=bar7\n"
        expected = expected1 + expected2 + expected3 + expected4 + expected5
        actual = test_cds2.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEquals(expected, actual)
        # what if name=None?
        test_cds2.name=[]
        expected1 = "sctg_0080_0020\tmaker\tCDS\t3734\t4034\t.\t+\t0\tID=foo1;Parent=bar7\n"
        expected2 = "sctg_0080_0020\tmaker\tCDS\t4092\t4332\t.\t+\t2\tID=foo2;Parent=bar7\n"
        expected3 = "sctg_0080_0020\tmaker\tCDS\t4399\t5185\t.\t+\t1\tID=foo3;Parent=bar7\n"
        expected4 = "sctg_0080_0020\tmaker\tCDS\t5249\t6565\t.\t+\t0\tID=foo4;Parent=bar7\n"
        expected5 = "sctg_0080_0020\tmaker\tCDS\t6630\t7436\t.\t+\t0\tID=foo5;Parent=bar7\n"
        expected = expected1 + expected2 + expected3 + expected4 + expected5
        actual = test_cds2.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEquals(expected, actual)

        # TODO test to_tbl when add annotations stuff

class TestExon(unittest.TestCase):

    def setUp(self):
        test_identifier1 = 3
        test_name1 = "BDOR_007864-RA:exon:0"
        test_indices1 = [3734, 4034] 
        test_score1 = 0.9
        test_parent_id1 = 2
        self.test_exon0 = Exon(identifier=test_identifier1, name=test_name1, indices=test_indices1, score=test_score1, parent_id=test_parent_id1)
        self.extra_identifiers = [4, 5, 6, 7]
        self.extra_names = ["BDOR_007864-RA:exon:1", "BDOR_007864-RA:exon:2", "BDOR_007864-RA:exon:3", "BDOR_007864-RA:exon:4"]
        self.extra_scores = [0.9, 0.9, 0.9, 0.9]
        self.extra_indices = [[4092, 4332], [4399, 5185], [5249, 6565], [6630, 7436]]
        self.test_exon1 = Exon(identifier=test_identifier1, name=test_name1, indices=test_indices1, score=test_score1, parent_id=test_parent_id1)
        for ind_pair in self.extra_indices:
            self.test_exon1.add_indices(ind_pair)
        for name in self.extra_names:
            self.test_exon1.add_name(name)
        for ident in self.extra_identifiers:
            self.test_exon1.add_identifier(ident)
        for score in self.extra_scores:
            self.test_exon1.add_score(score)


    def test_constructor(self):
        self.assertEquals('Exon', self.test_exon1.__class__.__name__)

    def test_add_indices(self):
        for ind_pair in self.extra_indices:
            self.test_exon0.add_indices(ind_pair)
        self.assertEquals(5, len(self.test_exon0.indices))
        self.assertEquals([5249, 6565], self.test_exon0.indices[3])

    def test_add_name(self):
        for name in self.extra_names:
            self.test_exon0.add_name(name)
        self.assertEquals(5, len(self.test_exon0.name))
        expected = "BDOR_007864-RA:exon:4"
        self.assertEquals(expected, self.test_exon0.name[4])

    def test_add_identifier(self):
        for ident in self.extra_identifiers:
            self.test_exon0.add_identifier(ident)
        self.assertEquals(5, len(self.test_exon0.identifier))
        self.assertEquals(6, self.test_exon0.identifier[3])

    def test_add_score(self):
        for score in self.extra_scores:
            self.test_exon0.add_score(score)
        self.assertEquals(5, len(self.test_exon0.score))
        self.assertEquals(0.9, self.test_exon0.score[4])
       
    def test_length(self): 
        self.assertEquals(301, self.test_exon0.length())
        self.assertEquals(3453, self.test_exon1.length())
       
    def test_adjust_indices(self): 
        self.test_exon1.adjust_indices(7)
        self.assertEquals(4339, self.test_exon1.indices[1][1])
        # adjust them back...
        self.test_exon1.adjust_indices(-7)
        self.assertEquals(4332, self.test_exon1.indices[1][1])

    def test_to_gff(self):
        expected1 = "sctg_0080_0020\tmaker\texon\t3734\t4034\t0.9\t+\t.\tID=3;Name=BDOR_007864-RA:exon:0;Parent=2\n"
        expected2 = "sctg_0080_0020\tmaker\texon\t4092\t4332\t0.9\t+\t.\tID=4;Name=BDOR_007864-RA:exon:1;Parent=2\n"
        expected3 = "sctg_0080_0020\tmaker\texon\t4399\t5185\t0.9\t+\t.\tID=5;Name=BDOR_007864-RA:exon:2;Parent=2\n"
        expected4 = "sctg_0080_0020\tmaker\texon\t5249\t6565\t0.9\t+\t.\tID=6;Name=BDOR_007864-RA:exon:3;Parent=2\n"
        expected5 = "sctg_0080_0020\tmaker\texon\t6630\t7436\t0.9\t+\t.\tID=7;Name=BDOR_007864-RA:exon:4;Parent=2\n"
        expected = expected1 + expected2 + expected3 + expected4 + expected5
        actual = self.test_exon1.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEquals(expected, actual)

        

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGenePart))
    suite.addTest(unittest.makeSuite(TestCDS))
    suite.addTest(unittest.makeSuite(TestExon))
    return suite

if __name__ == '__main__':
    unittest.main()

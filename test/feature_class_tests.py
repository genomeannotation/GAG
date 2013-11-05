#!/usr/bin/env python

import unittest
from mock import Mock
from src.feature_classes import GenePart, CDS, Exon, MRNA, Gene

class TestFeatureClasses(unittest.TestCase):

    def test_GenePart(self):
        # test constructor
        gp1 = GenePart()
        self.assertEqual(0, len(gp1.indices))
        self.assertFalse(gp1.score)
        self.assertFalse(gp1.parent_id)
        gp2 = GenePart(feature_type='CDS', indices=[1, 44])
        self.assertEqual(1, len(gp2.indices))
        self.assertEqual('CDS', gp2.feature_type)

        # test .add_indices
        gp2.add_indices([65, 103])
        self.assertEqual(2, len(gp2.indices))
        gp3 = GenePart(feature_type='exon')
        self.assertEqual(0, len(gp3.indices))
        gp3.add_indices([77, 144])
        self.assertEqual(1, len(gp3.indices))
        # error check
        self.assertRaises(ValueError, gp3.add_indices, 7)
        self.assertRaises(ValueError, gp3.add_indices, 'foo')

        # test .add_name
        self.assertEqual(0, len(gp3.name))
        gp3.add_name('BDOR_007864-RA:cds:0')
        self.assertEqual(1, len(gp3.name))

        # test .add_identifier
        self.assertEqual(0, len(gp3.identifier))
        gp3.add_identifier('7')
        self.assertEqual(1, len(gp3.identifier))

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
        gp2.identifier = ['foo1', 'foo2']
        gp2.parent_id = 'mama'
        expected = "ID=foo2;Parent=mama\n"
        self.assertEqual(expected, gp2.generate_attribute_entry(1))
        gp2.name = ['cds:0', 'cds:1']
        expected = "ID=foo1;Name=cds:0;Parent=mama\n"
        self.assertEqual(expected, gp2.generate_attribute_entry(0))
        # what if index out of range?
        self.assertFalse(gp2.generate_attribute_entry(2))
        # what if no identifier or parent_id?
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
        test_indices1 = [3734, 4034]
        extra_indices = [[4092, 4332], [4399, 5185], [5249, 6565], [6630, 7436]]
        test_name1 = "BDOR_007864-RA:cds:0"
        extra_names = ["BDOR_007864-RA:cds:1", "BDOR_007864-RA:cds:2", "BDOR_007864-RA:cds:3", "BDOR_007864-RA:cds:4"]
        test_identifier1 = 8
        extra_identifiers = [9, 10, 11, 12]
        test_phase1 = 0
        extra_phases = [2, 1, 0, 0]
        test_parent_id1 = 2
        test_cds1 = CDS(identifier=test_identifier1, name=test_name1, indices=test_indices1, score=None, phase=test_phase1, parent_id=test_parent_id1)
        self.assertEqual('CDS', test_cds1.__class__.__name__)
        # should also be able to construct w/o all the params...
        empty_cds = CDS()
        self.assertEqual('CDS', empty_cds.feature_type) 

        # test .add_indices
        for ind_pair in extra_indices:
            test_cds1.add_indices(ind_pair)
        self.assertEqual([4399, 5185], test_cds1.indices[2])

        # test .add_name
        for name in extra_names:
            test_cds1.add_name(name)
        self.assertEqual("BDOR_007864-RA:cds:4", test_cds1.name[4])

        # test .add_identifier
        for ident in extra_identifiers:
            test_cds1.add_identifier(ident)

        # test .add_phase
        for phase in extra_phases:
            test_cds1.add_phase(phase)

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
        # what if identifier, parent_id are strings? does it matter?
        test_cds2 = CDS(identifier='foo1', name=test_name1, indices=test_indices1, score=None, phase=test_phase1, parent_id='bar7')
        extra_identifiers2 = ['foo2', 'foo3', 'foo4', 'foo5']
        for ind_pair in extra_indices:
            test_cds2.add_indices(ind_pair)
        for name in extra_names:
            test_cds2.add_name(name)
        for ident in extra_identifiers2:
            test_cds2.add_identifier(ident)
        for phase in extra_phases:
            test_cds2.add_phase(phase)
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

        # TODO test to_tbl when add annotations stuff

    def test_Exon(self):
        # test constructor
        test_identifier1 = 3
        extra_identifiers = [4, 5, 6, 7]
        test_name1 = "BDOR_007864-RA:exon:0"
        extra_names = ["BDOR_007864-RA:exon:1", "BDOR_007864-RA:exon:2", "BDOR_007864-RA:exon:3", "BDOR_007864-RA:exon:4"]
        test_indices1 = [3734, 4034] 
        extra_indices = [[4092, 4332], [4399, 5185], [5249, 6565], [6630, 7436]]
        test_score1 = 0.9
        extra_scores = [0.9, 0.9, 0.9, 0.9]
        test_parent_id1 = 2
        test_exon1 = Exon(identifier=test_identifier1, name=test_name1, indices=test_indices1, score=test_score1, parent_id=test_parent_id1)
        self.assertEqual('Exon', test_exon1.__class__.__name__)

        # test .add_indices
        for ind_pair in extra_indices:
            test_exon1.add_indices(ind_pair)
        self.assertEqual(5, len(test_exon1.indices))

        # test .add_name
        for name in extra_names:
            test_exon1.add_name(name)

        # test .add_identifier
        for ident in extra_identifiers:
            test_exon1.add_identifier(ident)

        # test .add_score
        for score in extra_scores:
            test_exon1.add_score(score)
        
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
        test_mrna1 = MRNA(identifier=2, name="BDOR_007864-RA", indices=[3734, 7436], parent_id=1)
        self.assertEqual('MRNA', test_mrna1.__class__.__name__)

        # test .length
        self.assertEqual(3703, test_mrna1.length())

        # our mrna needs a cds and exon...
        fake_cds = Mock()
        fake_exon = Mock()

        # ... and maybe a start codon for fun
        fake_start_codon = Mock()

        # test .set_exon
        self.assertFalse(test_mrna1.exon)
        test_mrna1.set_exon(fake_exon)
        self.assertTrue(test_mrna1.exon)

        # test .set_cds
        self.assertFalse(test_mrna1.cds)
        test_mrna1.set_cds(fake_cds)
        self.assertTrue(test_mrna1.cds)
        
        # test .add_other_feature
        self.assertEquals(0, len(test_mrna1.other_features))
        test_mrna1.add_other_feature(fake_start_codon)
        self.assertEquals(1, len(test_mrna1.other_features))

        # test .adjust_indices
        test_mrna1.adjust_indices(32)
        self.assertEqual(7468, test_mrna1.indices[1])
        fake_cds.adjust_indices.assert_called_with(32)
        fake_exon.adjust_indices.assert_called_with(32)
        fake_start_codon.adjust_indices.assert_called_with(32)
        # adjust 'em back
        test_mrna1.adjust_indices(-32)
        self.assertEqual(3734, test_mrna1.indices[0])
        fake_cds.adjust_indices.assert_called_with(-32)
        fake_exon.adjust_indices.assert_called_with(-32)
        fake_start_codon.adjust_indices.assert_called_with(-32)

        # test .length_of_shortest_cds_segment
        fake_cds.length_of_shortest_segment.return_value = 241
        self.assertEquals(241, test_mrna1.length_of_shortest_cds_segment())
        fake_cds.length_of_shortest_segment.assert_called_once_with()

        # test .has_start, .has_stop
        fake_start_codon.feature_type = 'start_codon'
        self.assertTrue(test_mrna1.has_start())
        self.assertFalse(test_mrna1.has_stop())

        # test .to_gff (this is where the money is)
        fake_exon.to_gff.return_value = "...exon to gff\n"
        fake_cds.to_gff.return_value = "...cds to gff\n"
        fake_start_codon.to_gff.return_value = "...start codon to gff\n"
        expected = "sctg_0080_0020\tmaker\tmRNA\t"
        expected += "3734\t7436\t.\t+\t.\t"
        expected += "ID=2;Name=BDOR_007864-RA;Parent=1\n"
        expected += "...exon to gff\n...cds to gff\n"
        expected += "...start codon to gff\n"
        actual = test_mrna1.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEquals(expected, actual)
        fake_exon.to_gff.assert_called_with("sctg_0080_0020", "maker", '+')
        fake_cds.to_gff.assert_called_with("sctg_0080_0020", "maker", '+')
        fake_start_codon.to_gff.assert_called_with("sctg_0080_0020", "maker", '+')


    def test_Gene(self):
        # test constructor
        test_gene1 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[3734, 7436], strand='+', identifier=1, name="BDOR_007864")
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
        fake_mrna1.length_of_shortest_cds_segment.assert_called_with()
        fake_mrna2.length_of_shortest_cds_segment.assert_called_with()

        # test .adjust_indices
        test_gene1.adjust_indices(16)
        fake_mrna1.adjust_indices.assert_called_with(16)
        self.assertEquals(3750, test_gene1.indices[0])
        # adjust them back
        test_gene1.adjust_indices(-16)
        fake_mrna1.adjust_indices.assert_called_with(-16)
        self.assertEquals(3734, test_gene1.indices[0])
        # error check
        self.assertRaises(IndexError, test_gene1.adjust_indices, -4000)

        self.assertEquals(test_gene1.collidesRange(3500, 3800), True)
        self.assertEquals(test_gene1.collidesRange(3500, 3600), False)

        # test .to_gff
        fake_mrna1.to_gff.return_value = "fake mrna1 to gff here:)\n"
        fake_mrna2.to_gff.return_value = "fake mrna2 to gff here:)\n"
        expected = "sctg_0080_0020\tmaker\tgene\t3734\t7436\t.\t+\t."
        expected += "\tID=1;Name=BDOR_007864\n"
        expected += "fake mrna1 to gff here:)\n"
        expected += "fake mrna2 to gff here:)\n"
        self.assertEquals(expected, test_gene1.to_gff())

    def test_trim_end(self):
        # above tests are probably too long; should probably be
        # atomic. live and learn. this = trying to do better
        # test CDS.remove_segment
        cds = CDS(identifier="foo_CDS1", name="cds1", indices=[40, 75], phase=0, parent_id='mrna7')
        cds.add_identifier('foo_CDS2')
        cds.add_identifier('foo_CDS3')
        cds.add_name('cds2')
        cds.add_name('cds3')
        cds.add_indices([100, 175])
        cds.add_indices([200, 250])
        # CDS has 3 segments
        self.assertEquals(3, len(cds.indices))
        self.assertEquals(3, len(cds.name))
        self.assertEquals(3, len(cds.identifier))
        cds.remove_segment(2)
        # should now have 2 segments
        self.assertEquals(2, len(cds.indices))
        self.assertEquals(2, len(cds.name))
        self.assertEquals(2, len(cds.identifier))

        # test CDS.trim_end
        # (same idea for exon)
        cds.trim_end(150)
        # should keep both segments of cds
        self.assertEquals(2, len(cds.indices))
        self.assertEquals(2, len(cds.name))
        # and second segment's end index should be 230
        self.assertEquals(150, cds.indices[1][1])
        # now trim entire segment ...
        cds.trim_end(90)
        self.assertEquals(1, len(cds.indices))
        self.assertEquals(1, len(cds.name))
        # verify that all is okay if we trim off multiple segments
        cds.add_identifier('foo_CDS2')
        cds.add_identifier('foo_CDS3')
        cds.add_name('cds2')
        cds.add_name('cds3')
        cds.add_indices([100, 175])
        cds.add_indices([200, 250])
        self.assertEquals(3, len(cds.indices))
        self.assertEquals(3, len(cds.identifier))
        cds.trim_end(50)
        self.assertEquals(1, len(cds.indices))
        self.assertEquals(1, len(cds.name))
        self.assertEquals(50, cds.indices[0][1])
        self.assertEquals(40, cds.indices[0][0])

        # test GenePart.trim_end for one-segment feature
        stop_codon = GenePart(feature_type='stop_codon', identifier='foo_stop', name='stop1', indices=[20, 22], parent_id='foo_mrna')
        self.assertEquals(1, len(stop_codon.indices))
        self.assertEquals(22, stop_codon.indices[0][1])
        stop_codon.trim_end(21)
        self.assertEquals(1, len(stop_codon.indices))
        self.assertEquals(21, stop_codon.indices[0][1])
        stop_codon.trim_end(10)
        self.assertEquals(0, len(stop_codon.indices))
        self.assertEquals(0, len(stop_codon.identifier))

        # now do the same for mRNA
        mrna = MRNA(identifier='foo_mrna', name='mrna1', indices=[100, 200], parent_id='gene1')
        mrna.trim_end(150)
        self.assertEquals(150, mrna.indices[1])
        # now trim it into oblivion
        self.assertTrue(mrna.name)
        self.assertTrue(mrna.indices)
        mrna.trim_end(50)
        print(mrna.name)
        self.assertFalse(mrna.name)
        self.assertFalse(mrna.indices)
       
         # verify recursive call on child features
        mrna2 = MRNA(identifier='foo_mrna', name='mrna1', indices=[100, 200], parent_id='gene1')
        fake_exon = Mock()
        mrna2.set_exon(fake_exon)
        mrna2.trim_end(150)
        fake_exon.trim_end.assert_called_with(150)
        

        # ... and for Gene
        gene = Gene(seq_name="sctg_foo", source='maker', indices=[100, 200], strand='-', identifier='foo_gene', name='gene1')
        gene.trim_end(150)
        self.assertEquals(150, gene.indices[1])
        # now trim it into oblivion
        gene.trim_end(50)
        self.assertFalse(gene.indices)
        self.assertFalse(gene.identifier)
        
        # verify recursive call on child mrnas
        gene2 = Gene(seq_name="sctg_foo", source='maker', indices=[100, 200], strand='-', identifier='foo_gene', name='gene1')
        fake_mrna1 = Mock()
        fake_mrna2 = Mock()
        gene2.add_mrna(fake_mrna1)
        gene2.add_mrna(fake_mrna2)
        gene2.trim_end(150)
        fake_mrna1.trim_end.assert_called_with(150)
        fake_mrna2.trim_end.assert_called_with(150)

    def test_trim_begin(self):
        # test on exon; cds is about the same ...
        exon = Exon(identifier='foo_exon', name='exon1', indices=[50, 100], parent_id='mrna1')
        exon.add_identifier('foo_exon2')
        exon.add_name('exon2')
        exon.add_indices([200, 300])
        exon.add_identifier('foo_exon3')
        exon.add_name('exon3')
        exon.add_indices([400, 500])
        exon.trim_begin(20)
        self.assertEquals(31, exon.indices[0][0])
        self.assertEquals(381, exon.indices[2][0])
        # can result in negative indices; requires post-validation
        exon.trim_begin(50)
        self.assertEquals(-18, exon.indices[0][0])

        # verify same effect on one-segement other_features...
        utr = GenePart(feature_type='five_prime_utr', identifier='foo_utr', name='utr1', indices=[20, 22], parent_id='foo_mrna')
        utr.trim_begin(5)
        self.assertEquals(16, utr.indices[0][0])
        utr.trim_begin(20)
        self.assertEquals(-3, utr.indices[0][0])
        
        # do the same for mrna 
        mrna = MRNA(identifier='foo_mrna', name='mrna1', indices=[100, 200], parent_id='gene1')
        mrna.trim_begin(80)
        self.assertEquals(21, mrna.indices[0])
        # verify recursive call
        fake_cds = Mock()
        mrna.set_cds(fake_cds)
        mrna.trim_begin(10)
        self.assertEquals(12, mrna.indices[0])
        fake_cds.adjust_indices.assert_called_with(-9)
        # indices can reach negative territory
        mrna.trim_begin(50)
        self.assertEquals(-37, mrna.indices[0])
        
        # do the same for gene
        gene = Gene(seq_name="sctg_foo", source='maker', indices=[100, 200], strand='-', identifier='foo_gene', name='gene1')
        gene.trim_begin(80)
        self.assertEquals(21, gene.indices[0])
        # verify recursive call
        fake_mrna1 = Mock()
        fake_mrna2 = Mock()
        gene.add_mrna(fake_mrna1)
        gene.add_mrna(fake_mrna2)
        gene.trim_begin(10)
        self.assertEquals(112, gene.indices[1])
        # note -- gene doesn't call trim_begin on child features;
        # it calls adjust_indices. which means that GenePart.trim_begin
        # is probably useless. live and learn.
        fake_mrna1.adjust_indices.assert_called_with(-9)
        fake_mrna2.adjust_indices.assert_called_with(-9)


        
        
        
        

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFeatureClasses))
    return suite

if __name__ == '__main__':
    unittest.main()

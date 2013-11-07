#!/usr/bin/env python

import unittest
from mock import Mock, PropertyMock
from src.feature_classes import GenePart, CDS, Exon, MRNA, Gene

class TestMRNA(unittest.TestCase):

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
        cds.trim_end(150)
        self.assertEquals(0, cds.indices[2][0])
        self.assertEquals(0, cds.indices[2][1])
        self.assertEquals(150, cds.indices[1][1])
        cds.trim_end(90)
        self.assertEquals(0, cds.indices[1][0])
        self.assertEquals(0, cds.indices[1][1])
        self.assertEquals(75, cds.indices[0][1]) #unchanged

        # test GenePart.trim_end for one-segment feature
        stop_codon = GenePart(feature_type='stop_codon', identifier='foo_stop', name='stop1', indices=[20, 22], parent_id='foo_mrna')
        self.assertEquals(22, stop_codon.indices[0][1])
        stop_codon.trim_end(21)
        self.assertEquals(21, stop_codon.indices[0][1])
        stop_codon.trim_end(10)
        self.assertEquals(0, stop_codon.indices[0][0])
        self.assertEquals(0, stop_codon.indices[0][0])

        # now do the same for mRNA
        mrna = MRNA(identifier='foo_mrna', name='mrna1', indices=[100, 200], parent_id='gene1')
        mrna.trim_end(150)
        self.assertEquals(100, mrna.indices[0])
        self.assertEquals(150, mrna.indices[1])
        mrna.trim_end(50)
        self.assertEquals(0, mrna.indices[0])
        self.assertEquals(0, mrna.indices[1])
       
         # verify recursive call on child features
        mrna2 = MRNA(identifier='foo_mrna', name='mrna1', indices=[100, 200], parent_id='gene1')
        fake_exon = Mock()
        mrna2.set_exon(fake_exon)
        mrna2.trim_end(150)
        fake_exon.trim_end.assert_called_with(150)
        

        # ... and for Gene
        gene = Gene(seq_name="sctg_foo", source='maker', indices=[100, 200], strand='-', identifier='foo_gene', name='gene1')
        gene.trim_end(150)
        self.assertEquals(100, gene.indices[0])
        self.assertEquals(150, gene.indices[1])
        gene.trim_end(50)
        self.assertEquals(0, gene.indices[0])
        self.assertEquals(0, gene.indices[1])
        
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
        # only need to test on gene since gene will call
        # adjust_indices recursively (this is already under test)
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
        fake_mrna1.adjust_indices.assert_called_with(-9)
        fake_mrna2.adjust_indices.assert_called_with(-9)

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

        # test on MRNA
        mrna = MRNA(identifier='foo', name='foo', indices=[1, 100], parent_id='bar')
        fake_cds = Mock()
        mrna.set_cds(fake_cds)
        mrna.adjust_phase()
        fake_cds.adjust_phase.assert_called_with()

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

        # test on mRNA...
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

        # test on Gene
        nice_gene = Gene(seq_name='fooseq', source='maker', indices=[-23, 127], strand='-', identifier='foo', name='foo')
        junk_gene = Gene(seq_name='barseq', source='maker', indices=[-400, -100], strand='-', identifier='bar', name='bar')
        nice_gene.clean_up_indices()
        self.assertEquals(1, nice_gene.indices[0])
        self.assertEquals(127, nice_gene.indices[1])
        junk_gene.clean_up_indices()
        self.assertEquals(0, junk_gene.indices[0])
        self.assertEquals(0, junk_gene.indices[1])
        # test recursive call
        mrna1 = Mock()
        mrna2 = Mock()
        nice_gene.add_mrna(mrna1)
        nice_gene.add_mrna(mrna2)
        nice_gene.clean_up_indices()
        mrna1.clean_up_indices.assert_called_with()
        mrna2.clean_up_indices.assert_called_with()
        
    def test_remove_trimmed_segments(self):
        exon = Exon(identifier='foo', name='foo', indices=[0, 0], parent_id='foo')
        exon.add_identifier('foo2')
        exon.add_name('foo2')
        exon.add_indices([22, 100])
        self.assertEquals(2, len(exon.identifier))
        self.assertEquals(2, len(exon.indices))
        exon.remove_trimmed_segments()
        self.assertEquals(1, len(exon.name))
        self.assertEquals(1, len(exon.indices))
        # verify that it will remove all segments if appropriate
        exon.indices[0] = [0, 0]
        exon.remove_trimmed_segments()
        self.assertEquals(0, len(exon.indices))
        self.assertEquals(0, len(exon.identifier))

    def test_valid_codon(self):
        stop = GenePart(feature_type='stop_codon', identifier='foo', name='foo', indices=[2, 3], parent_id='foofoo')
        self.assertFalse(stop.valid_codon())

    def test_remove_invalid_features(self):
        mrna = MRNA(identifier='foo', name='foo', indices=[20, 50], parent_id='foo') 
        nice_exon = Mock()
        empty_cds = Mock() 
        empty_stop_codon = Mock()
        invalid_start_codon = Mock()
        nice_utr = Mock()
        mrna.set_exon(nice_exon)
        mrna.set_cds(empty_cds)
        mrna.add_other_feature(empty_stop_codon)
        mrna.add_other_feature(nice_utr)
        mrna.add_other_feature(invalid_start_codon)
        nice1 = PropertyMock(return_value = [10, 50])
        nice2 = PropertyMock(return_value = [1, 9])
        invalid1 = PropertyMock(return_value = [1, 2])
        empty = PropertyMock(return_value=[])
        starttype = PropertyMock(return_value = 'start_codon')
        invalid_start_codon.valid_codon.return_value = False
        type(nice_exon).indices = nice1
        type(nice_utr).indices = nice2
        type(invalid_start_codon).indices = invalid1
        type(invalid_start_codon).feature_type = starttype
        type(empty_cds).indices = empty
        type(empty_stop_codon).indices = empty
        self.assertEquals(3, len(mrna.other_features))
        self.assertTrue(mrna.cds)
        self.assertTrue(mrna.exon)
        mrna.remove_invalid_features()
        self.assertFalse(mrna.cds)
        self.assertEquals(1, len(mrna.other_features))

        # try it on Gene
        gene = Gene(seq_name="sctg_foo", source='maker', indices=[100, 200], strand='-', identifier='foo_gene', name='gene1')
        junk_mrna = MRNA(identifier='junk', name='junk', indices=[0, 0], parent_id='foo_gene')
        nice_mrna = Mock()
        nice_ind = PropertyMock(return_value = [100, 200])
        type(nice_mrna).indices = nice_ind
        gene.add_mrna(junk_mrna)
        gene.add_mrna(nice_mrna) 
        self.assertEquals(2, len(gene.mrnas))
        gene.remove_invalid_features()
        self.assertEquals(1, len(gene.mrnas))
        nice_mrna.remove_invalid_features.assert_called_with()

    def test_trim(self):
        gene = Gene(seq_name='sctg_foo', source='maker', indices=[100, 200], strand='-', identifier='foo_gene', name='gene1')
        nice_mrna = Mock()
        bad_mrna = Mock()
        nice_ind = PropertyMock(return_value = [100, 200])
        bad_ind = PropertyMock(return_value = [0, 0])
        type(nice_mrna).indices = nice_ind
        type(bad_mrna).indices = bad_ind
        gene.add_mrna(nice_mrna)
        gene.add_mrna(bad_mrna)
        bed_indices = [120, 180]
        self.assertEquals(2, len(gene.mrnas))
        gene.trim(bed_indices)
        self.assertEquals(1, len(gene.mrnas))
        self.assertEquals(1, gene.indices[0])
        self.assertEquals(61, gene.indices[1])
        nice_mrna.trim_end.assert_called_with(180)
        bad_mrna.trim_end.assert_called_with(180)
        nice_mrna.adjust_indices.assert_called_with(-119)
        bad_mrna.adjust_indices.assert_called_with(-119)
        nice_mrna.adjust_phase.assert_called_with()
        bad_mrna.adjust_phase.assert_called_with()
        nice_mrna.clean_up_indices.assert_called_with()
        bad_mrna.clean_up_indices.assert_called_with()
        nice_mrna.remove_invalid_features.assert_called_with()

        
        
        
        

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMRNA))
    return suite

if __name__ == '__main__':
    unittest.main()

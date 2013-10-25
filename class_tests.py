#!/usr/bin/env python

import unittest
from feature_classes import CDS

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
        result = test_cds1.to_gff(seq_name="sctg_0080_0020", source="maker", strand='+')
        self.assertEqual(expected, result)

##########################
if __name__ == '__main__':
    unittest.main()

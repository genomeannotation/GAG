#!/usr/bin/env python

import unittest
from src.stats_manager import StatsManager

class TestStatsManager(unittest.TestCase):

    def setUp(self):
        self.mgr = StatsManager()

    def test_initialize(self):
        self.assertEquals(self.mgr.ref_stats["num_CDS"], 0)

    def test_clear_ref(self):
        self.populate_ref()
        self.assertEquals(self.mgr.ref_stats["num_CDS"], 7)
        self.mgr.clear_ref()
        self.assertEquals(self.mgr.ref_stats["num_CDS"], 0)


    def populate_ref(self):
        self.mgr.ref_stats["seq_length"] = 100
        self.mgr.ref_stats["num_genes"] = 5
        self.mgr.ref_stats["num_mRNA"] = 7
        self.mgr.ref_stats["num_CDS"] = 7
        self.mgr.ref_stats["longest_gene"] = 25
        self.mgr.ref_stats["longest_mRNA"] = 25
        self.mgr.ref_stats["longest_CDS"] = 20
        self.mgr.ref_stats["shortest_gene"] = 10
        self.mgr.ref_stats["shortest_mRNA"] = 10
        self.mgr.ref_stats["shortest_CDS"] = 6
        self.mgr.ref_stats["total_gene_length"] = 70
        self.mgr.ref_stats["total_mRNA_length"] = 70
        self.mgr.ref_stats["total_CDS_length"] = 60

    def get_new_dict(self):
        d = {}
        d["seq_length"] = 50
        d["num_genes"] = 1
        d["num_mRNA"] = 1
        d["num_CDS"] = 1
        d["longest_gene"] = 30
        d["longest_mRNA"] = 30
        d["longest_CDS"] = 8
        d["shortest_gene"] = 5
        d["shortest_mRNA"] = 5
        d["shortest_CDS"] = 3
        d["total_gene_length"] = 15
        d["total_mRNA_length"] = 15
        d["total_CDS_length"] = 10
        return d

        
    def test_update_ref(self):
        self.populate_ref()
        newdict = self.get_new_dict()
        self.assertEquals(self.mgr.ref_stats["seq_length"], 100)
        self.assertEquals(self.mgr.ref_stats["shortest_CDS"], 6)
        self.assertEquals(self.mgr.ref_stats["longest_gene"], 25)
        self.mgr.update_ref(newdict)
        self.assertEquals(self.mgr.ref_stats["seq_length"], 150)
        self.assertEquals(self.mgr.ref_stats["shortest_CDS"], 3)
        self.assertEquals(self.mgr.ref_stats["longest_gene"], 30)

    def test_summary(self):
        self.populate_ref()
        expected = "\t\tReference Genome\tModified Genome\n"
        expected += "\t\t----------------\t---------------\n"
        expected += "seq_length:\t\t100\t\t0\n"
        expected += "num_genes:\t\t5\t\t0\n"
        expected += "num_mRNA:\t\t7\t\t0\n"
        expected += "num_CDS:\t\t7\t\t0\n"
        # The next three lines have longer stat-names,
        # so fewer tabs so things will line up.
        # Is there a smarter way to do that?
        expected += "total_gene_length:\t70\t\t0\n"
        expected += "total_mRNA_length:\t70\t\t0\n"
        expected += "total_CDS_length:\t60\t\t0\n"
        expected += "shortest_gene:\t\t10\t\t0\n"
        expected += "shortest_mRNA:\t\t10\t\t0\n"
        expected += "shortest_CDS:\t\t6\t\t0\n"
        expected += "longest_gene:\t\t25\t\t0\n"
        expected += "longest_mRNA:\t\t25\t\t0\n"
        expected += "longest_CDS:\t\t20\t\t0\n"
        summary = self.mgr.summary()
        self.assertEquals(summary, expected)


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestStatsManager))
    return suite

if __name__ == '__main__':
    unittest.main()

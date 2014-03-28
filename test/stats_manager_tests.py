#!/usr/bin/env python

import unittest
from src.stats_manager import StatsManager

class TestStatsManager(unittest.TestCase):

    def setUp(self):
        self.mgr = StatsManager()

    def test_initialize(self):
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


''' 
    def stats(self):
        stats = dict()
        
        stats["seq_length"] = len(self.bases)
        stats["num_genes"] = len(self.genes)
        stats["num_mRNA"] = self.get_num_mrna()
        stats["num_CDS"] = self.get_num_cds()
        stats["longest_gene"] = self.get_longest_gene()
        stats["longest_mRNA"] = self.get_longest_mrna()
        stats["longest_CDS"] = self.get_longest_cds()
        stats["shortest_gene"] = self.get_shortest_gene()
        stats["shortest_mRNA"] = self.get_shortest_mrna()
        stats["shortest_CDS"] = self.get_shortest_cds()
        stats["total_mRNA_length"] = self.get_total_mrna_length()
        stats["total_CDS_length"] = self.get_total_cds_length()
'''


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestStatsManager))
    return suite

if __name__ == '__main__':
    unittest.main()

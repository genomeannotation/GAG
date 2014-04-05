#!/usr/bin/env python

import unittest
from src.stats_manager import StatsManager
from src.stats_manager import format_column
from src.stats_manager import format_columns

class TestStatsManager(unittest.TestCase):

    def setUp(self):
        self.mgr = StatsManager()

    def test_initialize(self):
        self.assertEquals(self.mgr.ref_stats["num_CDS"], 0)

    def test_clear_alt(self):
        self.mgr.update_alt(self.get_new_dict())
        self.assertEquals(self.mgr.alt_stats["num_CDS"], 1)
        self.mgr.clear_alt()
        self.assertEquals(self.mgr.alt_stats["num_CDS"], 0)

    def populate_ref(self):
        self.mgr.ref_stats["seq_length"] = 100
        self.mgr.ref_stats["num_genes"] = 5
        self.mgr.ref_stats["num_mRNA"] = 7
        self.mgr.ref_stats["num_exons"] = 7
        self.mgr.ref_stats["num_introns"] = 7
        self.mgr.ref_stats["num_CDS"] = 7
        self.mgr.ref_stats["CDS: complete"] = 3
        self.mgr.ref_stats["CDS: start, no stop"] = 1
        self.mgr.ref_stats["CDS: stop, no start"] = 1
        self.mgr.ref_stats["CDS: no stop, no start"] = 2
        self.mgr.ref_stats["longest_gene"] = 25
        self.mgr.ref_stats["longest_mRNA"] = 25
        self.mgr.ref_stats["longest_exon"] = 21
        self.mgr.ref_stats["longest_intron"] = 21
        self.mgr.ref_stats["longest_CDS"] = 20
        self.mgr.ref_stats["shortest_gene"] = 10
        self.mgr.ref_stats["shortest_mRNA"] = 10
        self.mgr.ref_stats["shortest_exon"] = 8
        self.mgr.ref_stats["shortest_intron"] = 8
        self.mgr.ref_stats["shortest_CDS"] = 6
        self.mgr.ref_stats["total_gene_length"] = 70
        self.mgr.ref_stats["total_mRNA_length"] = 70
        self.mgr.ref_stats["total_exon_length"] = 65
        self.mgr.ref_stats["total_intron_length"] = 65
        self.mgr.ref_stats["total_CDS_length"] = 60

    def get_new_dict(self):
        d = {}
        d["seq_length"] = 50
        d["num_genes"] = 1
        d["num_mRNA"] = 1
        d["num_exons"] = 1
        d["num_introns"] = 1
        d["num_CDS"] = 1
        d["CDS: complete"] = 3
        d["CDS: start, no stop"] = 1
        d["CDS: stop, no start"] = 1
        d["CDS: no stop, no start"] = 2
        d["longest_gene"] = 30
        d["longest_mRNA"] = 30
        d["longest_exon"] = 9
        d["longest_intron"] = 9
        d["longest_CDS"] = 8
        d["shortest_gene"] = 5
        d["shortest_mRNA"] = 5
        d["shortest_exon"] = 2
        d["shortest_intron"] = 2
        d["shortest_CDS"] = 3
        d["total_gene_length"] = 15
        d["total_mRNA_length"] = 15
        d["total_exon_length"] = 15
        d["total_intron_length"] = 15
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
        expected =  "                                     Reference Genome     Modified Genome     \n"
        expected += "                                     ----------------     ---------------     \n"
        expected += "seq_length                           100                  0                   \n"
        expected += "num_genes                            5                    0                   \n"
        expected += "num_mRNA                             7                    0                   \n"
        expected += "num_exons                            7                    0                   \n"
        expected += "num_introns                          7                    0                   \n"
        expected += "num_CDS                              7                    0                   \n"
        expected += "CDS: complete                        3                    0                   \n"
        expected += "CDS: start, no stop                  1                    0                   \n"
        expected += "CDS: stop, no start                  1                    0                   \n"
        expected += "CDS: no stop, no start               2                    0                   \n"
        expected += "total_gene_length                    70                   0                   \n"
        expected += "total_mRNA_length                    70                   0                   \n"
        expected += "total_exon_length                    65                   0                   \n"
        expected += "total_intron_length                  65                   0                   \n"
        expected += "total_CDS_length                     60                   0                   \n"
        expected += "shortest_gene                        10                   0                   \n"
        expected += "shortest_mRNA                        10                   0                   \n"
        expected += "shortest_exon                        8                    0                   \n"
        expected += "shortest_intron                      8                    0                   \n"
        expected += "shortest_CDS                         6                    0                   \n"
        expected += "longest_gene                         25                   0                   \n"
        expected += "longest_mRNA                         25                   0                   \n"
        expected += "longest_exon                         21                   0                   \n"
        expected += "longest_intron                       21                   0                   \n"
        expected += "longest_CDS                          20                   0                   \n"
        expected += "mean gene length                     14.0                 0                   \n"
        expected += "mean mRNA length                     10.0                 0                   \n"
        expected += "mean exon length                     9.28571428571        0                   \n"
        expected += "mean intron length                   9.28571428571        0                   \n"
        expected += "mean CDS length                      8.57142857143        0                   \n"
        expected += "prop. of genome covered by genes     0.7                  0                   \n"
        expected += "prop. of genome covered by CDS       0.6                  0                   \n"
        expected += "mRNAs per gene                       1.4                  0                   \n"
        expected += "exons per mRNA                       1.0                  0                   \n"
        expected += "introns per mRNA                     1.0                  0                   \n"
        summary = self.mgr.summary()
        self.assertEquals(summary, expected)
        
    def test_format_column(self):
        column = ['a', 'sd', 'asdf']
        self.assertEquals(format_column(column, 5), ['a        ', 'sd       ', 'asdf     '])
        
    def test_format_columns(self):
        desired_tbl = '    columnA columnB \n' \
                      '    ------- ------- \n' \
                      'dog 24      4222    \n' \
                      'foo 4232234 84      \n'
        column_names = ['columnA', 'columnB']
        dictA = {'foo' : 4232234, 'dog' : 24}
        dictB = {'foo' : 84, 'dog' : 4222}
        self.assertEquals(format_columns(column_names, ['dog', 'foo'], [dictA, dictB], 1), desired_tbl)


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestStatsManager))
    return suite

if __name__ == '__main__':
    unittest.main()

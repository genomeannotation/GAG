#!/usr/bin/env python
# coding=utf-8

import unittest
from src.stats_manager import StatsManager
from src.stats_manager import format_column
from src.stats_manager import format_columns


class TestStatsManager(unittest.TestCase):
    def setUp(self):
        self.mgr = StatsManager()

    def test_initialize(self):
        self.assertEquals(self.mgr.ref_stats["Number of CDS"], 0)

    def test_clear_alt(self):
        self.mgr.update_alt(self.get_new_dict())
        self.assertEquals(self.mgr.alt_stats["Number of CDS"], 1)
        self.mgr.clear_alt()
        self.assertEquals(self.mgr.alt_stats["Number of CDS"], 0)

    def test_clear_all(self):
        self.populate_ref()
        self.mgr.update_alt(self.get_new_dict())
        self.assertEquals(self.mgr.alt_stats["Number of CDS"], 1)
        self.assertEquals(self.mgr.ref_stats["Number of CDS"], 7)
        self.mgr.clear_all()
        self.assertEquals(self.mgr.alt_stats["Number of CDS"], 0)
        self.assertEquals(self.mgr.ref_stats["Number of CDS"], 0)

    def populate_ref(self):
        self.mgr.ref_stats["Total sequence length"] = 100
        self.mgr.ref_stats["Number of genes"] = 5
        self.mgr.ref_stats["Number of mRNAs"] = 7
        self.mgr.ref_stats["Number of exons"] = 7
        self.mgr.ref_stats["Number of introns"] = 7
        self.mgr.ref_stats["Number of CDS"] = 7
        self.mgr.ref_stats["Overlapping genes"] = 3
        self.mgr.ref_stats["Contained genes"] = 3
        self.mgr.ref_stats["CDS: complete"] = 3
        self.mgr.ref_stats["CDS: start, no stop"] = 1
        self.mgr.ref_stats["CDS: stop, no start"] = 1
        self.mgr.ref_stats["CDS: no stop, no start"] = 2
        self.mgr.ref_stats["Longest gene"] = 25
        self.mgr.ref_stats["Longest mRNA"] = 25
        self.mgr.ref_stats["Longest exon"] = 21
        self.mgr.ref_stats["Longest intron"] = 21
        self.mgr.ref_stats["Longest CDS"] = 20
        self.mgr.ref_stats["Shortest gene"] = 10
        self.mgr.ref_stats["Shortest mRNA"] = 10
        self.mgr.ref_stats["Shortest exon"] = 8
        self.mgr.ref_stats["Shortest intron"] = 8
        self.mgr.ref_stats["Shortest CDS"] = 6
        self.mgr.ref_stats["Total gene length"] = 70
        self.mgr.ref_stats["Total mRNA length"] = 70
        self.mgr.ref_stats["Total exon length"] = 65
        self.mgr.ref_stats["Total intron length"] = 65
        self.mgr.ref_stats["Total CDS length"] = 60

    @staticmethod
    def get_new_dict():
        d = {"Total sequence length": 50,
             "Number of genes": 1,
             "Number of mRNAs": 1,
             "Number of exons": 1,
             "Number of introns": 1,
             "Number of CDS": 1,
             "Overlapping genes": 1,
             "Contained genes": 1,
             "CDS: complete": 3,
             "CDS: start, no stop": 1,
             "CDS: stop, no start": 1,
             "CDS: no stop, no start": 2,
             "Longest gene": 30,
             "Longest mRNA": 30,
             "Longest exon": 9,
             "Longest intron": 9,
             "Longest CDS": 8,
             "Shortest gene": 5,
             "Shortest mRNA": 5,
             "Shortest exon": 2,
             "Shortest intron": 2,
             "Shortest CDS": 3,
             "Total gene length": 15,
             "Total mRNA length": 15,
             "Total exon length": 15,
             "Total intron length": 15,
             "Total CDS length": 10}
        return d

    def test_alt_is_empty(self):
        self.assertTrue(self.mgr.alt_is_empty())
        self.mgr.update_alt(self.get_new_dict())
        self.assertFalse(self.mgr.alt_is_empty())

    def test_update_ref(self):
        self.populate_ref()
        newdict = self.get_new_dict()
        self.assertEquals(self.mgr.ref_stats["Total sequence length"], 100)
        self.assertEquals(self.mgr.ref_stats["Shortest CDS"], 6)
        self.assertEquals(self.mgr.ref_stats["Longest gene"], 25)
        self.mgr.update_ref(newdict)
        self.assertEquals(self.mgr.ref_stats["Total sequence length"], 150)
        self.assertEquals(self.mgr.ref_stats["Shortest CDS"], 3)
        self.assertEquals(self.mgr.ref_stats["Longest gene"], 30)

    def test_summary_with_modifications(self):
        self.populate_ref()
        self.mgr.update_alt(self.get_new_dict())
        expected = "                                 Reference Genome     Modified Genome     \n"
        expected += "                                 ----------------     ---------------     \n"
        expected += "Total sequence length            100                  50                  \n"
        expected += "Number of genes                  5                    1                   \n"
        expected += "Number of mRNAs                  7                    1                   \n"
        expected += "Number of exons                  7                    1                   \n"
        expected += "Number of introns                7                    1                   \n"
        expected += "Number of CDS                    7                    1                   \n"
        expected += "Overlapping genes                3                    1                   \n"
        expected += "Contained genes                  3                    1                   \n"
        expected += "CDS: complete                    3                    3                   \n"
        expected += "CDS: start, no stop              1                    1                   \n"
        expected += "CDS: stop, no start              1                    1                   \n"
        expected += "CDS: no stop, no start           2                    2                   \n"
        expected += "Total gene length                70                   15                  \n"
        expected += "Total mRNA length                70                   15                  \n"
        expected += "Total exon length                65                   15                  \n"
        expected += "Total intron length              65                   15                  \n"
        expected += "Total CDS length                 60                   10                  \n"
        expected += "Shortest gene                    10                   5                   \n"
        expected += "Shortest mRNA                    10                   5                   \n"
        expected += "Shortest exon                    8                    2                   \n"
        expected += "Shortest intron                  8                    2                   \n"
        expected += "Shortest CDS                     6                    3                   \n"
        expected += "Longest gene                     25                   30                  \n"
        expected += "Longest mRNA                     25                   30                  \n"
        expected += "Longest exon                     21                   9                   \n"
        expected += "Longest intron                   21                   9                   \n"
        expected += "Longest CDS                      20                   8                   \n"
        expected += "mean gene length                 14                   15                  \n"
        expected += "mean mRNA length                 10                   15                  \n"
        expected += "mean exon length                 9                    15                  \n"
        expected += "mean intron length               9                    15                  \n"
        expected += "mean CDS length                  9                    10                  \n"
        expected += "% of genome covered by genes     70.0                 30.0                \n"
        expected += "% of genome covered by CDS       60.0                 20.0                \n"
        expected += "mean mRNAs per gene              1                    1                   \n"
        expected += "mean exons per mRNA              1                    1                   \n"
        expected += "mean introns per mRNA            1                    1                   \n"
        summary = self.mgr.summary()
        self.assertEquals(summary, expected)

    def test_summary_without_modifications(self):
        self.mgr.update_alt(self.get_new_dict())
        self.mgr.update_ref(self.get_new_dict())
        expected = "                                 Genome     \n" + \
                   "                                 ------     \n" + \
                   "Total sequence length            50         \n" + \
                   "Number of genes                  1          \n" + \
                   "Number of mRNAs                  1          \n" + \
                   "Number of exons                  1          \n" + \
                   "Number of introns                1          \n" + \
                   "Number of CDS                    1          \n" + \
                   "Overlapping genes                1          \n" + \
                   "Contained genes                  1          \n" + \
                   "CDS: complete                    3          \n" + \
                   "CDS: start, no stop              1          \n" + \
                   "CDS: stop, no start              1          \n" + \
                   "CDS: no stop, no start           2          \n" + \
                   "Total gene length                15         \n" + \
                   "Total mRNA length                15         \n" + \
                   "Total exon length                15         \n" + \
                   "Total intron length              15         \n" + \
                   "Total CDS length                 10         \n" + \
                   "Shortest gene                    5          \n" + \
                   "Shortest mRNA                    5          \n" + \
                   "Shortest exon                    2          \n" + \
                   "Shortest intron                  2          \n" + \
                   "Shortest CDS                     3          \n" + \
                   "Longest gene                     30         \n" + \
                   "Longest mRNA                     30         \n" + \
                   "Longest exon                     9          \n" + \
                   "Longest intron                   9          \n" + \
                   "Longest CDS                      8          \n" + \
                   "mean gene length                 15         \n" + \
                   "mean mRNA length                 15         \n" + \
                   "mean exon length                 15         \n" + \
                   "mean intron length               15         \n" + \
                   "mean CDS length                  10         \n" + \
                   "% of genome covered by genes     30.0       \n" + \
                   "% of genome covered by CDS       20.0       \n" + \
                   "mean mRNAs per gene              1          \n" + \
                   "mean exons per mRNA              1          \n" + \
                   "mean introns per mRNA            1          \n"
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
        dict_a = {'foo': 4232234, 'dog': 24}
        dict_b = {'foo': 84, 'dog': 4222}
        self.assertEquals(format_columns(column_names, ['dog', 'foo'], [dict_a, dict_b], 1), desired_tbl)


def suite():
    _suite = unittest.TestSuite()
    _suite.addTest(unittest.makeSuite(TestStatsManager))
    return _suite


if __name__ == '__main__':
    unittest.main()

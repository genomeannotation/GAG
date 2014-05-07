#!/usr/bin/env python

import unittest
from mock import Mock
from src.filter_manager import FilterManager

class TestFilterManager(unittest.TestCase):

    def setUp(self):
        self.filter_mgr = FilterManager()
        
    def test_modify_filter_arg(self):
        self.filter_mgr.set_filter_arg('cds_shorter_than', '30')
        self.assertEqual(self.filter_mgr.get_filter_arg('cds_shorter_than'), 30)




##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFilterManager))
    return suite

if __name__ == '__main__':
    unittest.main()

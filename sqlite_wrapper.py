#!/usr/bin/env python

import sqlite3

class SqliteWrapper:

	def __init__(self, db):
		self.database = db

	def getRowAsStr(self, table, id_qual, row_id):
		c = self.database.cursor()
		c.execute('SELECT * FROM '+table+' WHERE '+id_qual+'='+row_id+' LIMIT 1')
		return str(c.fetchone())

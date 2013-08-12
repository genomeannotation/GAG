#!/usr/bin/env python

import sqlite3

class SqliteWrapper:

	def __init__(self, path):
		self.database = sqlite3.connect(path)

	def createTable(self, table, columns):
		c = self.database.cursor()
		c.execute('CREATE TABLE '+table+'('+columns+')')

	def insertRow(self, table, columns):
		c = self.database.cursor()
		c.execute('INSERT INTO '+table+' VALUES('+columns+')')

	def commit():
		self.database.commit()

	def getRowAsStr(self, table, id_qual, row_id):
		c = self.database.cursor()
		c.execute('SELECT * FROM '+table+' WHERE '+id_qual+'='+row_id+' LIMIT 1')
		return str(c.fetchone())

#!/usr/bin/env python

import sqlite3

class SqliteWrapper:

	def __init__(self, db):
		self.database = db

	def createTable(self, table, columns):
		cols = ''
		for val in columns:
			cols = cols+str(val)+','
		cols=cols[:len(cols)-1]

		c = self.database.cursor()
		c.execute('CREATE TABLE '+table+'('+cols+')')

	def insertRow(self, table, columns):
		cols = ''
		for val in columns:
			cols = cols+str(val)+','
		cols=cols[:len(cols)-1]

		c = self.database.cursor()
		c.execute('INSERT INTO '+table+' VALUES('+cols+')')

	def commit():
		self.database.commit()

	def getRowAsStr(self, table, id_qual, row_id):
		c = self.database.cursor()
		c.execute('SELECT * FROM '+table+' WHERE '+id_qual+'='+row_id+' LIMIT 1')
		return str(c.fetchone())

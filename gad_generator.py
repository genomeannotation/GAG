#!/usr/bin/env python

import sqlite3

print('Hello, world!')

# Connect to the database.
db = sqlite3.connect('test_files/example.db')

# Grab the database's cursor to execute commands
c = db.cursor()

# Create the table
c.execute('CREATE TABLE Cars(Id integer primary key, Name text, Cost integer)')
cur.execute("INSERT INTO Cars VALUES(1,'Audi',52642)")

# Commit the changes to the database
db.commit()

# Close the database when we're done.
db.close()

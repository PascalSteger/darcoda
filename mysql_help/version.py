#!/usr/bin/python
# -*- coding: utf-8 -*-

import MySQLdb
import sys
import initialize as my

connection = MySQLdb.connect('igloo.dhcp.phys','psteger','Pinux10','astro')
cursor = connection.cursor ()

cursor.execute ("SELECT VERSION()")

row = cursor.fetchone ()
print 'Server version:', row[0]

cursor.close ()
connection.close ()

sys.exit()

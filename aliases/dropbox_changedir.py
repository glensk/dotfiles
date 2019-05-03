#!/usr/bin/env python

# Script for showing or setting the dropbox folder.
#
# Execute without options to show current dropbox folder (if non-default).
# Execute with --setfolder=/foo/bar to set new dropbox folder.
#
# I dedicate this work to the public domain by waiving all my rights to the
# work worldwide under copyright law, including all related and neighboring
# rights, to the extent allowed by law.
#
# You can copy, modify, distribute and perform the work, even for commercial
# purposes, all without asking permission.
#
# Wim Coenen (wcoenen@gmail.com).

import base64
import optparse
import os
import os.path
import sqlite3

# parse command line options
cmdparser = optparse.OptionParser()
cmdparser.add_option("-s","--setfolder", dest="folder",
  help="set dropbox folder")
(options, args) = cmdparser.parse_args()

db_path = os.path.expanduser("~/.dropbox/dropbox.db")
db = sqlite3.connect(db_path)
cursor = db.cursor()

# get dropbox_path
cursor.execute("select value from config where key='dropbox_path'")
dropbox_path = "<default>"
for entry in cursor:
   dropbox_path_base64 = entry[0]
   dropbox_path_raw = base64.decodestring(dropbox_path_base64)
   dropbox_path = dropbox_path_raw[1:-5]
print "current dropbox path: %s" % dropbox_path

if not options.folder is None:
   new_path_raw = "V" + os.path.abspath(options.folder) + "\np1\n."
   new_path_base64 = base64.encodestring(new_path_raw) 
   cursor.execute("delete from config where key='dropbox_path'")
   cursor.execute("insert into config (key,value) values (?,?)", \
      ("dropbox_path", new_path_base64))
   db.commit()
   print "new dropbox path: %s" % options.folder
   
db.close()

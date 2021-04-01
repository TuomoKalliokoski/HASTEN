# 
# Copyright (c) 2021 Orion Corporation
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimer in the documentation 
# and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software 
# without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Simulate docking i.e. read score from a textfile and insert empty pose to
# database
#
# param 1 "confs". here just compound IDs
# param 2 database name
# param 3 dock.txt (the complete docking results)
# param 4 hasten<->smilesid mapping
# param 5 iteration 
#
import sys
import csv
import sqlite3

iteration = int(sys.argv[5])
alldocks = {}
# read in docking scores from text file
with open(sys.argv[3]) as dockfile:
    for row in csv.reader(dockfile,delimiter=" "):
        alldocks[row[1]] = float(row[0])

# read in codes and hasten ids
smilesids_to_hastenids = {}
with open(sys.argv[4]) as hastenfile:
    for row in csv.reader(hastenfile,delimiter=" "):
        r = row[0].split("|")
        if r[0] not in smilesids_to_hastenids:
            smilesids_to_hastenids[r[0]] = r[1]
try:
    conn=sqlite3.connect(sys.argv[2])
except:
    print("Error while accessing database!")
    sys.exit(1)
finally:
    if conn:
        c=conn.cursor()
        to_db = []
        for smilesid in smilesids_to_hastenids: 
            to_db.append((alldocks[smilesid],smilesids_to_hastenids[smilesid]))
        c.executemany("UPDATE data SET dock_score = ? WHERE hastenid = ?",to_db)
        conn.commit()
        to_db = []
        for smilesid in smilesids_to_hastenids: 
            to_db.append((iteration,smilesids_to_hastenids[smilesid]))
        c.executemany("UPDATE data SET dock_iteration = ? WHERE hastenid = ?",to_db)
        conn.commit()
        conn.close()

        print("simulate_docking.py OK")

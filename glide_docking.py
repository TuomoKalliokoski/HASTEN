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
# parse Glide output for HASTEN
#
# 1 = docked poses in maestro format
# 2 = confs to be docked in maestro format
# 3 = failed dock score value
# 4 = database name
# 5 = smilesid to hastenid mapping
# 6 = cutoff for storing posees
#
# Written with Schrodinger Suite 2020-3
#
import sys
import csv
import sqlite3
from schrodinger import structure

# read in hastenids (in docked file they use smilesids)
smilesid_to_hastenid = {}
for line in open(sys.argv[5],"rt"):
    l = line.strip().split("|")
    smilesid_to_hastenid[l[0]] = l[1]

score_cutoff = float(sys.argv[6])

mols = {}
poses = {}
reader = structure.StructureReader(sys.argv[1])
for st in reader:
    s=st.property["s_m_title"]
    hastenid=smilesid_to_hastenid[s]
    docking_score=float(st.property["r_i_docking_score"])
    if hastenid not in mols or mols[hastenid]>=docking_score:
        mols[hastenid] = docking_score
        if docking_score <= score_cutoff:
            poses[hastenid] = "{ \n  s_m_m2io_version\n  :::\n  2.0.0 \n} \n\n"+structure.write_ct_to_string(st)

# check not docked mols and give them bad score
reader = structure.StructureReader(sys.argv[2])
bad_score = float(sys.argv[3])
for st in reader:
    s=st.property["s_m_title"]
    hastenid=smilesid_to_hastenid[s]
    if hastenid not in mols:
        mols[hastenid] = bad_score

try:
    conn=sqlite3.connect(sys.argv[4])
except:
    print("Error while accesing database")
    sys.exit(1)
finally:
    if conn:
        c=conn.cursor()
        to_db = []
        for hastenid in mols:
            to_db.append((mols[hastenid],hastenid))
        c.executemany("UPDATE data SET dock_score = ? WHERE hastenid = ?",to_db)
        conn.commit()
        to_db = []
        for hastenid in poses:
            to_db.append((hastenid,bytes(poses[hastenid],encoding="utf-8")))
        c.executemany("REPLACE INTO poses(hastenid,pose) VALUES (?,?)",to_db)
        conn.commit()
        conn.close()
        print("glide_docking.py done.")

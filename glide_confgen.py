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
# parse LigPrep output for HASTEN
#
# parameters:
# 1 = .maegz output from LigPrep
# 2 = .db HASTEN database
#
# Written with Schrodinger Suite 2020-3
#
import sys
import csv
import sqlite3
from schrodinger import structure

mols = {}
reader = structure.StructureReader(sys.argv[1])
for st in reader:
    s=st.property["s_m_title"].split("|")
    st.property["s_m_title"]=s[0]
    if s[1] not in mols:
        mols[s[1]] = "{ \n  s_m_m2io_version\n  :::\n  2.0.0 \n} \n\n"
    mols[s[1]]+=structure.write_ct_to_string(st)

try:
    conn=sqlite3.connect(sys.argv[2])
except:
    print("Error while accessing database!")
finally:
    if conn:
        c=conn.cursor()
        to_db = []
        for hastenid in mols:
            to_db.append((hastenid,bytes(mols[hastenid],encoding="utf-8")))
        c.executemany("REPLACE INTO confs(hastenid,conf) VALUES (?,?)",to_db)
        conn.commit()
        conn.close()
        print("glide_confgen.py done.")

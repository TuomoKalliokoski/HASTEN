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
# "Generate" conformations for database
#
# Find out first internal hastenid for your molecules and then
# add correct binary blob to confs table (dummy structure with just
# molecule name)
#
#
import sys
import csv
import sqlite3

mols = []
confs = []
with open(sys.argv[1]) as smilesfile:
    for row in csv.reader(smilesfile,delimiter=" "):
        mols.append(row[1].split("|")[1])
        conf=bytes(row[1].split("|")[0]+"\n",encoding="utf-8")
        confs.append(conf)
try:
    conn=sqlite3.connect(sys.argv[2])
except:
    print("Error while accessing database")
    sys.exit(1)
finally:
    if conn:
        c=conn.cursor()
        # insert empty mols for these
        to_db = []
        for counter in range(len(mols)): 
            to_db.append((mols[counter],confs[counter]))
        c.executemany("REPLACE INTO confs(hastenid,conf) VALUES (?,?)",to_db)
        conn.commit()
        conn.close()
        print("simulate_confgen.py done.")

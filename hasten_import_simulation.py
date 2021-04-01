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
"""
HASTEN simulation import

Import testing data for HASTEN
"""


import argparse
import os
import sys
import csv
import sqlite3

def parse_cmd_line():
    """
    Parse command line using ArgumentParser

    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(description="Import data for HASTEN simulation")
    parser.add_argument("-s","--smiles",required=True,type=str,help="SMILES input file")
    parser.add_argument("-d","--dock",required=True,type=str,help="Docking score file")
    parser.add_argument("-o","--output",required=True,type=str,help="Output database filename")
    return parser.parse_args()

def files_exist(args):
    """
    Check if input files exist

    :param args: parsed arguments
    :return: True if files exist, false otherwise
    """
    if not os.path.exists(args.smiles):
        print("SMILES file missing")
        return False
    if not os.path.exists(args.dock):
        print("Docking score file missing")
        return False
    if os.path.exists(args.output):
        os.remove(args.output)
    return True

def check_files(mols,docks):
    """
    Verify that we have docking score for each SMILES

    :param mols:Dictionary of SMILES
    :param docks:Dictionar of scores
    :return:True if the data is OK, false otherwise
    """
    for molid in mols:
        if molid not in docks:
            return False
    return True

def import_db(args):
    """
    Process the files and import them to database

    :param args: Parsed arguments
    """
    mols = {}
    docks = {}
    with open(args.smiles) as smilesfile:
        for row in csv.reader(smilesfile,delimiter=" "):
            mols[row[1]] = row[0]
    with open(args.dock) as dockfile:
        for row in csv.reader(dockfile,delimiter=" "):
            try:
                docks[row[1]] = float(row[0])
            except:
                print("Invalid numeric value in docking scores:",row)
                sys.exit(1)
    if len(mols)!=len(docks) or not check_files(mols,docks):
        print("Error: SMILES and docking scores do not match")
        sys.exit(1)
    try:
        conn=sqlite3.connect(args.output)
    except error in e:
        print("Error while creating database",e)
    finally:
        if conn:
            c=conn.cursor()
            c.execute("CREATE TABLE data (hastenid INTEGER PRIMARY KEY,smiles TEXT,smilesid TEXT,dock_score NUMERIC,dock_iteration INTEGER,pred_score NUMERIC,dataset_status INTEGER)")
            c.execute("CREATE TABLE confs (hastenid INTEGER PRIMARY KEY,conf BLOB)")
            c.execute("CREATE TABLE poses (hastenid INTEGER PRIMARY KEY,pose BLOB)")
            to_db = []
            for molname in mols: to_db.append((mols[molname],molname))
            c.executemany("INSERT INTO data(smiles,smilesid) VALUES (?,?)",to_db)
            conn.commit()
            conn.close()
    
if __name__ == "__main__":
    args = parse_cmd_line()
    if (not files_exist(args)): sys.exit(1)
    import_db(args)

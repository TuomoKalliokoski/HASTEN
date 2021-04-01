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
HASTEN analyze simulation

You can calculate recalls using this if you have the whole docking data for
your datbase (like in the benchmarking study)
"""
import argparse
import os
import sys
import csv
import sqlite3
import subprocess

def parse_cmd_line():
    """
    Parse command line using ArgumentParser

    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(description="Analyze HASTEN simulation")
    parser.add_argument("-m","--database",required=True,type=str,help="HASTEN database")
    parser.add_argument("-d","--dock",required=True,type=str,help="Docking data")
    return parser.parse_args()

def files_exist(args):
    """
    Check if input files exist

    :param args: parsed arguments
    :return: True if files exist, false otherwise
    """
    if not os.path.exists(args.database):
        print("HASTEN database file missing!")
        return False
    if not os.path.exists(args.dock):
        print("Docking data missing!")
        return False
    return True

def analyze(args):
    """
    Calculate recalls

    :param args: Parsed arguments
    """
    # read in docking scores in ugly way as the dock.txt can be huge
    # assuming it is sorted! 
    docks=int(subprocess.check_output("/usr/bin/wc -l "+args.dock,shell=True).split()[0])
    one_percent=round(0.01*docks)
    row_count = 0
    with open(args.dock) as dockfile:
        for row in csv.reader(dockfile,delimiter=" "):
            try:
                if row_count>=one_percent:
                    cutoff = float(row[0])
                    break
                row_count += 1
            except:
                print("Invalid numeric value in docking scores")
                sys.exit(1)
    print("WARNING: Assuming that docking scores are sorted in the file.")
    print()
    print("DOCKING DATA")
    print("Docked compounds:",docks)
    print("Top 1% of docked compounds:",one_percent)
    print("Docking score cutoff at top 1%:",cutoff)
    print()
    print("HASTEN DATA\n")
    print("Iter\t Mols\t Recall")
    print("-------------------------------------")
    
    # sort by docking score
    # get the top 1% docking score cut off and number of hits
    sqlstr="SELECT dock_score,dock_iteration FROM data WHERE data.dock_score<="+str(cutoff)
    try:
        conn=sqlite3.connect(args.database)
    except error in e:
        print("Error while accessing database:",e)
        sys.exit(1)
    if conn:
        hasten_data = {}
        c = conn.cursor()
        cur=c.execute(sqlstr)
        for row in cur.fetchmany(4000000):
            if row[1] not in hasten_data:
                hasten_data[row[1]] = 0
            hasten_data[row[1]] += 1
        conn.close()
        so_far = 0
        for iteration in sorted(hasten_data.keys()):
            so_far += hasten_data[iteration]
            print(iteration,"\t",hasten_data[iteration],"\t",round(so_far/one_percent,3))

if __name__ == "__main__":
    args = parse_cmd_line()
    if (not files_exist(args)): sys.exit(1)
    analyze(args)

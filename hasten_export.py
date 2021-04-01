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
HASTEN export

Export results from HASTEN runs
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
    parser = argparse.ArgumentParser(description="Export data from HASTEN")
    parser.add_argument("-m","--database",required=True,type=str,help="HASTEN database")
    parser.add_argument("-c","--cutoff",required=True,type=float,help="Docking score cut off for molecules")
    parser.add_argument("-z","--out-dock-poses",required=False,type=str,help="Docked poses from docking calculations")
    parser.add_argument("-x","--out-dock-scores",required=False,type=str,help="Scores for the docked poses")
    parser.add_argument("-a","--out-pred-confs",required=False,type=str,help="Conformers for to be docked compounds")
    parser.add_argument("-q","--out-pred-scores",required=False,type=str,help="Predicted scores for the to be docked poses")
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
    return True

def export_to_file(args,outputtype):
    """
    Output data to file

    :param args: parsed arguments
    :param outputtype: "dock-poses","dock-scores","pred-confs","pred-scores"
    """
    if outputtype=="dock-poses":
        outputfile=args.out_dock_poses
        sqlstr="SELECT pose FROM poses INNER JOIN data ON data.hastenid==poses.hastenid WHERE data.dock_score<="+str(args.cutoff)
        output_blob = True
    elif outputtype=="dock-scores":
        outputfile=args.out_dock_scores
        sqlstr="SELECT smiles,smilesid,dock_score,dock_iteration FROM data WHERE data.dock_score<="+str(args.cutoff)
        output_blob = False
    elif outputtype=="pred-confs":
        outputfile=args.out_pred_confs
        sqlstr="SELECT conf FROM confs INNER JOIN data ON data.hastenid==confs.hastenid WHERE data.dock_score IS NULL AND data.pred_score<="+str(args.cutoff)
        output_blob = True
    elif outputtype=="pred-scores":
        outputfile=args.out_pred_scores
        sqlstr="SELECT smiles,smilesid,pred_score FROM data WHERE data.dock_score IS NULL AND data.pred_score<="+str(args.cutoff)
        output_blob = False
    else:
        print("INTERNAL ERROR: invalid export_to_file definition")
        sys.exit(1)
    print("Exporting to",outputfile)
    try:
        conn=sqlite3.connect(args.database)
    except error in e:
        print("Error while accessing database:",e)
        sys.exit(1)
    if conn:
        c = conn.cursor()
        if output_blob:
            rowblob=c.execute(sqlstr).fetchall()
            w = open(outputfile,"wb")
            for row in rowblob: w.write(row[0])
            w.close()
        else:
            rowtxt=c.execute(sqlstr).fetchall()
            with open(outputfile, "wt") as csvfile:
                reswriter=csv.writer(csvfile,delimiter=",")
                for row in rowtxt: reswriter.writerow(row)
        conn.close()

if __name__ == "__main__":
    args = parse_cmd_line()
    if (not files_exist(args)): sys.exit(1)
    if args.out_dock_poses is not None: export_to_file(args,"dock-poses")
    if args.out_dock_scores is not None: export_to_file(args,"dock-scores")
    if args.out_pred_confs is not None: export_to_file(args,"pred-confs")
    if args.out_pred_scores is not None: export_to_file(args,"pred-scores")


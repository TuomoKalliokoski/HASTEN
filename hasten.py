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
HASTEN main program

This takes care of most stuff outside database import.
"""

import argparse
import os
import sys
import csv
import sqlite3
import random
import tempfile
import glob

def parse_cmd_line():
    """
    Parse command line using ArgumentParser

    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(description="Run HASTEN")
    parser.add_argument("-m","--database",required=False,type=str,help="HASTEN database")
    parser.add_argument("-p","--protocol",required=True,type=str,help="Screening protocol file")
    parser.add_argument("-i","--iteration",required=False,type=int,help="Iteration number to start")

    parser.add_argument("-a","--hand-operate",required=False,type=str,choices=["dock","train","split-pred","pred","import-pred"],help="Hand-operated mode (only for expert users)")
    return parser.parse_args()

def files_exist(args):
    """
    Check if input files exist

    :param args: parsed arguments
    :return: True if files exist, false otherwise
    """
    if args.database is not None and not os.path.exists(args.database):
        print("HASTEN database file missing!")
        return False
    if not os.path.exists(args.protocol):
        print("Screening protocol file missing!")
        return False
    return True

def get_protocol(filename):
    """
    Load a protocol file

    :param filename: Filename for the .protocol file
    :return: Protocol dictionary
    """
    protocol = {}
    keywords = {"confgen":"file","docking":"file","ml_train":"file","ml_pred":"file","dataset_size":"float","dataset_split":"floats3","train_mode":"text","random_seed":"integer","pred_size":"integer","stop_criteria":"integer","pred_split":"integer"}
    for line in open(filename,"rt"):
        if line[0].startswith("#"):
            continue
        l = line.strip().split("=")
        if len(l)!=2 or l[0] not in keywords:
            print("Invalid line in protocol:"+line.strip())
            sys.exit(1)
        if keywords[l[0]] == "file":
            if not os.path.exists(l[1]):
                print("checking",l[0])
                print(l[0],"is defined but is missing in the protocol file")
                sys.exit(1)
            else:
                protocol[l[0]] = l[1]
        elif keywords[l[0]] == "float":
                try:
                    protocol[l[0]] = float(l[1])
                except:
                    print(l[0],"must be float in the protocol file")
                    sys.exit(1)
        elif keywords[l[0]] == "floats3":
            protocol[l[0]] = []
            numbers = l[1].split(" ")
            if len(numbers)!=3:
                print(l[0],"must be three floats seperated by space in the protocol file")
                sys.exit(1)
            for number in numbers:
                if len(number.split("."))!=2:
                    print(l[0],"must be three floats seperated by space in the protocol file")
                    sys.exit(1)
                else:
                    protocol[l[0]].append(float(number))
        elif keywords[l[0]] == "text":
            protocol[l[0]] = l[1]
        elif keywords[l[0]] == "integer":
            try:
                protocol[l[0]] = int(l[1])
            except:
                print(l[0],"needs to be an integer in the protocol file.")
                sys.exit(1)

    for keyword in keywords:
        if keyword not in protocol:
            print(keyword,"not defined in protocol.")
            sys.exit(1)

    return protocol

def pick_compounds_for_docking(protocol,db,iteration):
    """
    Pick set of compounds from database for docking and save them into
    output files.

    :param protocol: Protocol dictionary
    :param db: The filename of SQlite3 database
    :param iteration: iteration integer
    """
    try:
        conn=sqlite3.connect(db)
    except:
        print("Error while accessing database!")
        sys.exit(1)
    if conn:
        c = conn.cursor()
        number_of_mols=c.execute("SELECT COUNT() FROM data").fetchone()[0]
        number_to_dock=int(round(protocol["dataset_size"]*number_of_mols))
        print("Number of molecules in the database",number_of_mols)
        print("Picking",protocol["dataset_size"]*100,"% for docking (",number_to_dock,")...")

        # first time pick just random set (slow way but as this is done once in
        # every iteration it does not matter)
        if iteration==1:
            to_dock=c.execute("SELECT hastenid FROM data WHERE dock_score IS NULL ORDER BY random() LIMIT ?",[number_to_dock,]).fetchall()
        else:
            to_dock=c.execute("SELECT hastenid FROM data WHERE dock_score IS NULL ORDER BY pred_score LIMIT ?",[number_to_dock,]).fetchall()
        smilesids = []
        for smile_id in to_dock: smilesids.append(smile_id[0])
       
        # check by joining to confs table which of the mols have already confs 
        sqlstr="SELECT data.hastenid FROM data INNER JOIN confs ON confs.hastenid=data.hastenid WHERE confs.conf IS NOT NULL AND data.hastenid IN {}".format(str(tuple(smilesids))).replace("[","").replace("]","")
        confs_avail=c.execute(sqlstr).fetchall()
        conn.close()
   
        # add to list those that do not have conf 
        confgen_smilesids = []
        for smile_id in confs_avail: 
            confgen_smilesids.append(smile_id[0])
        calc_confs = []
        for smile_id in smilesids:
            if smile_id not in confgen_smilesids:
                calc_confs.append(smile_id)
    
        return(smilesids,calc_confs)

def run_confgen(protocol,db,smilesids):
    """
    Run outside conformer generator (simply starts external code)

    :param protocol: Protocol dictionary
    :param db: The filename of SQlite3 database
    :param smilesids: List of hastenids to be exported
    """
    try:
        conn=sqlite3.connect(db)
    except:
        print("Error while accessing database!")
        sys.exit(1)
    if conn:
        c = conn.cursor()
        sqlstr="SELECT smiles,smilesid,hastenid FROM data WHERE hastenid IN {}".format(str(tuple(smilesids))).replace("[","").replace("]","")
        rowsmiles=c.execute(sqlstr).fetchall()
        conn.close()
 
        temp_name = tempfile.mkstemp(".smi","hasten_confgen_","/tmp")[1]
        w = open(temp_name,"wt")
        for row in rowsmiles:
            w.write(row[0]+" "+row[1]+"|"+str(row[2])+"\n")
        w.close()
        os.system(protocol["confgen"]+" "+temp_name+" "+db)
        # clean up if the confgen script did not already
        try:
            os.unlink(temp_name)
        except:
            pass

def run_docking(protocol,db,smilesids,iteration):
    """
    Run outside docking (simply starts external code)

    :param protocol: Protocol dictionary
    :param db: The filename of SQlite3 database
    :param smilesids: List of hastenids to be exported
    :param iteration: iteration integer
    """
    db_chunk_size=123456
    try:
        conn=sqlite3.connect(db)
    except:
        print("Error while accessing database!")
        sys.exit(1)
    if conn:
        c = conn.cursor()
        sqlstr="SELECT conf,data.hastenid,data.smilesid FROM confs INNER JOIN data ON data.hastenid=confs.hastenid WHERE data.hastenid IN {}".format(str(tuple(smilesids))).replace("[","").replace("]","")
 
        temp_name = tempfile.mkstemp(".out","hasten_dock_confs_","/tmp")[1]
        temp2_name = tempfile.mkstemp(".txt","hasten_dock_ids_","/tmp")[1]
        w = open(temp_name,"wb")
        w2 = open(temp2_name,"wt")
        db_cursor = c.execute(sqlstr)
        rowsmiles=db_cursor.fetchmany(db_chunk_size)
        while(len(rowsmiles)>0):
            for row in rowsmiles: 
                w.write(row[0])
                w2.write(str(row[2])+"|"+str(row[1])+"\n")
            rowsmiles=db_cursor.fetchmany(db_chunk_size)
        w.close()
        w2.close()

        os.system(protocol["docking"]+" "+temp_name+" "+db+" "+temp2_name+" "+str(iteration))
        # clean up if the confgen script did not already
        conn.close()
        try:
            os.unlink(temp_name)
        except:
            pass
        try:
            os.unlink(temp2_name)
        except:
            pass

def run_ml_train(protocol,db,iteration):
    """
    Pick set of compounds from database for docking and save them into
    output files.

    :param protocol: Protocol dictionary
    :param db: The filename of SQlite3 database
    :param iteration: iteration integer
    """
    try:
        conn=sqlite3.connect(db)
    except error:
        print("Error while accessing database!")
        sys.exit(1)
    if conn:
        c = conn.cursor()
        sqlstr="SELECT smiles,hastenid,dock_score FROM data WHERE dock_score IS NOT NULL"
        rowsmiles=c.execute(sqlstr).fetchall()
        dataset_size=len(rowsmiles)
        print(dataset_size,"compounds with docking result")
        conn.close()
        if protocol["train_mode"]=="scratch":
            print("Runninng in scratch mode")
            # calculate number of molecules into each set
            train_set_size = int(round(protocol["dataset_split"][0]*dataset_size))
            validation_set_size = int(round(protocol["dataset_split"][1]*dataset_size))
            test_set_size = int(round(protocol["dataset_split"][2]*dataset_size))
            print("Training set:",train_set_size)
            print("Validation set:",validation_set_size)
            print("Test set:",test_set_size)
            # in scratch mode we don't care where everything goes
            random.shuffle(rowsmiles)
            train_set = []
            for i in range(train_set_size):
                train_set.append(rowsmiles.pop())
            validation_set = []
            for i in range(validation_set_size):
                validation_set.append(rowsmiles.pop())
            # dump the rest to testset
            test_set = []
            while(len(rowsmiles)>0): test_set.append(rowsmiles.pop())
        elif protocol["train_mode"] == "increase":
            print("FIXME: increase training mode not yet implemented!")
            sys.exit(1)
        else:
            print("Error, invalid train_mode in the protocol file:",protocol["train_mode"])
            sys.exit(1)
    train_filename = write_for_ml(train_set)
    valid_filename = write_for_ml(validation_set)
    test_filename = write_for_ml(test_set)

    os.system(protocol["ml_train"]+" "+train_filename+" "+valid_filename+" "+test_filename+" iter"+str(iteration))

    os.unlink(train_filename)
    os.unlink(valid_filename)
    os.unlink(test_filename)

def run_ml_pred(protocol,db,iteration,mode="normal"):
    """
    Predict compounds either in automatic or hand-operated mode (see mode)

    :param protocol: Protocol dictionary
    :param db: The filename of SQlite3 database
    :param iteration: iteration integer
    :param mode: string "split" means hand-operated split, "para" means prediction in hand-operated mode and default "normal" the automatic mode
    """
    if mode=="split":
        conn=sqlite3.connect(db)
        c = conn.cursor()
        sqlstr="SELECT COUNT(*) FROM data WHERE dock_score IS NULL"
        number_of_comps=int(c.execute(sqlstr).fetchall()[0][0])
        sqlstr="SELECT smiles,hastenid FROM data WHERE dock_score IS NULL"
        per_machine = int(number_of_comps/protocol["pred_split"])
        print("Compounds per machine optimal:",per_machine)
        cur_machine=1
        cur_machine_count=0
        cur_chunk = 1
        # only predict those that we don't have docking_score yet
        sqlstr="SELECT smiles,hastenid FROM data WHERE dock_score IS NULL"
        db_batch_size = 4000000
        db_cursor = c.execute(sqlstr)
        rowsmiles = c.fetchmany(db_batch_size)
        chunk = []
        while len(rowsmiles) > 0:
            while len(rowsmiles) > 0:
                chunk.append(rowsmiles.pop())
                cur_machine_count+=1
                if len(chunk)>=protocol["pred_size"] or cur_machine_count>=per_machine:
                    if not os.path.exists("PRED"+str(cur_machine)):
                        os.mkdir("PRED"+str(cur_machine))
                    chunk_filename = write_for_ml(chunk,with_score=False,filename="PRED"+str(cur_machine)+"/iter"+str(iteration)+"_pred_input_"+str(cur_machine)+"_"+str(cur_chunk)+".csv")
                    chunk = []
                    cur_chunk += 1
                    print("Wrote",chunk_filename)
                    if cur_machine_count>=per_machine:
                        cur_machine_count = 0
                        cur_machine += 1
                        cur_chunk = 1
                        if not os.path.exists("PRED"+str(cur_machine)):
                            os.mkdir("PRED"+str(cur_machine))
            rowsmiles = db_cursor.fetchmany(db_batch_size)
        conn.close()
        # write leftover compounds in last chunk
        if len(chunk)>0:
            chunk_filename = write_for_ml(chunk,with_score=False,filename="PRED"+str(cur_machine)+"/iter"+str(iteration)+"_pred_input_"+str(cur_machine)+"_"+str(cur_chunk+1)+".csv")
    elif mode=="para":
        for filename in glob.glob("iter*_pred_input_*.csv"):
            print("Predicting:",filename)
            pred_chunk(protocol,None,None,iteration,filename)
    elif mode=="normal":
        print("Running prediction in normal, single-CPU mode...")
        conn=sqlite3.connect(db)
        db_cursor = conn.cursor()
        sqlstr="SELECT smiles,hastenid FROM data WHERE dock_score IS NULL"
        # note we must load everything here to RAM as we are also updating
        # the DB inside the loop - sqlite3 gets stuck otherwise
        db_cursor = db_cursor.execute(sqlstr)
        rowsmiles = db_cursor.fetchall()
        conn.close()
        while len(rowsmiles)>0:
            chunk = []
            while len(chunk)<=protocol["pred_size"] and len(rowsmiles)>0:
                chunk.append(rowsmiles.pop())
            if len(chunk)>0:
                print(len(rowsmiles),"compounds left to be ranked by the ML model")
                pred_chunk(protocol,db,chunk,iteration)
    else:
        print("BUG IN ml_chemprop_pred(): invalid mode!")
        sys.exit(2)
        
def pred_chunk(protocol,db,chunk,iteration,filename=None):
    """
    Predict a chunk of molecules

    :param protocol: Protocol dictionary
    :param db: The filename of SQlite3 database
    :param chunk: List of hastenids to be predicted
    :param iteration: iteration integer
    :param filename: automatic mode is None, hand-operated mode has the filename
    """
    if filename is None:
        chunk_filename = write_for_ml(chunk,with_score=False)
        chunk_output = tempfile.mkstemp(".csv","hasten","/tmp")[1]
    else:
        chunk_filename = filename
        chunk_output = filename.replace("_input_","_output_")
    os.system(protocol["ml_pred"]+" "+chunk_filename+" iter"+str(iteration)+" "+chunk_output)
    if filename is None:
        write_pred_to_db(db,chunk_output)
        os.unlink(chunk_filename)
        os.unlink(chunk_output)

def write_pred_to_db(db,filename):
    """
    Write predictions to db from a ML output file

    :param db: The filename of SQlite3 database
    :param filename: The filename of the ML output file
    """
    pred_scores = []
    with open(filename) as outputfile:
        # skip header row
        csvreader = csv.reader(outputfile,delimiter=",")
        next(csvreader)
        for row in csvreader:
            pred_scores.append((float(row[2]),row[1]))
    try:
        conn=sqlite3.connect(db)
    except:
        print("Error while accesing database")
        sys.exit(1)
    finally:
        if conn:
            print("Importing predictions to db...")
            c=conn.cursor()
            c.executemany("UPDATE data SET pred_score = ? WHERE hastenid = ?",pred_scores)
            conn.commit()
            conn.close()

# with_score = do we have score or not
def write_for_ml(rows,with_score=True,filename=None):
    """
    Write data for ml training

    :param rows: The data rows to be written
    :param with_score: Write data with score (usually True)
    :param filename: If None, write into temporary file
    :return: Filename for the temporary file
    """
    if filename is None:
        temp_name = tempfile.mkstemp(".smi","hasten","/tmp")[1]
    else:
        temp_name = filename
    w = open(temp_name,"wt")
    if with_score:
        w.write("smiles,hastenid,docking_score\n")
    else:
        w.write("smiles,hastenid\n")
    for row in rows:
        if with_score:
            w.write(row[0]+","+str(row[1])+","+str(row[2])+"\n")
        else:
            w.write(row[0]+","+str(row[1])+"\n")
    w.close()
    return temp_name

def run_ml_import(protocol,db):
    """
    Import bunch of files in hand-operated mode from ML predictions

    :param protocol: Protocol dictionary
    :param db: The database filename
    """

    conn=sqlite3.connect(db)
    if conn:
        for filename in glob.glob("iter*_output_*.csv"):
            print("Importing predictions from",filename)
            write_pred_to_db(db,filename)

def run_hasten(protocol,args):
    """
    Starts the whole process

    :param protocol: Protocol dictionary
    :param args: Parsed arguments
    """
    random.seed(protocol["random_seed"])
    print()
    print("    )        (                 )  ")
    print(" ( /(  (     )\ )  *   )    ( /(  ")
    print(" )\()) )\   (()/(` )  /((   )\()) ")
    print("((_)((((_)(  /(_))( )(_))\ ((_)\  ")
    print(" _((_)\ _ )\(_)) (_(_()|(_) _((_) ")
    print("| || (_)_\(_) __||_   _| __| \| | ")
    print("| __ |/ _ \ \__ \  | | | _|| .` | ")
    print("|_||_/_/ \_\|___/  |_| |___|_|\_| ")
    print()
    print("HASTEN (macHine leArning booSTEd dockiNg) version 0.1")
    print("https://github.com/TuomoKalliokoski/HASTEN")
    if args.iteration is not None:
        iteration = args.iteration
    else:
        iteration = 1

    if args.hand_operate is not None:
        print("Hand-operated mode activated.")
        print("Iteration",iteration)
        if args.hand_operate == "dock":
            compounds_for_docking,compounds_for_confgen = pick_compounds_for_docking(protocol,args.database,iteration)
            print(len(compounds_for_confgen),"molecules to conformer generation...")
            run_confgen(protocol,args.database,compounds_for_confgen)
            print("Running docking...")
            run_docking(protocol,args.database,compounds_for_docking,iteration)
        elif args.hand_operate == "train":
            print("Running machine learning training...")
            run_ml_train(protocol,args.database,iteration)
        elif args.hand_operate == "split-pred":
            print("Splitting for machine learning predictions...")
            run_ml_pred(protocol,args.database,iteration,mode="split")
        elif args.hand_operate == "pred":
            print("Running machine learning predictions...")
            run_ml_pred(protocol,args.database,iteration,mode="para")
        elif args.hand_operate == "import-pred":
            print("Importing machine learning predictions...")
            run_ml_import(protocol,args.database)
    else:
        while iteration<=protocol["stop_criteria"]:
                print("Iteration",iteration)

                if iteration>1:
                    print("Running machine learning training...")
                    run_ml_train(protocol,args.database,iteration)
                    print("Running machine learning prediction...")
                    run_ml_pred(protocol,args.database,iteration)

                compounds_for_docking,compounds_for_confgen = pick_compounds_for_docking(protocol,args.database,iteration)

                print(len(compounds_for_confgen),"molecules to conformer generation...")
                run_confgen(protocol,args.database,compounds_for_confgen)
                print("Running docking...")
                run_docking(protocol,args.database,compounds_for_docking,iteration)

                iteration+=1

    print("\nHASTEN finished.")

if __name__ == "__main__":
    args = parse_cmd_line()
    if (not files_exist(args)): sys.exit(1)
    protocol=get_protocol(args.protocol)
    run_hasten(protocol,args)

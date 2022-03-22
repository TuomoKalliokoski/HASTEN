# 
# Copyright (c) 2021-2022 Orion Corporation
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

"""

import argparse
import os
import sys
import csv
import sqlite3
import random
import tempfile
import glob
import shutil
import statistics

def parse_cmd_line():
    """
    Parse command line using ArgumentParser

    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(description="Run HASTEN")
    parser.add_argument("-m","--database",required=False,type=str,help="HASTEN database")
    parser.add_argument("-i","--iteration",required=False,type=int,help="Iteration number")

    parser.add_argument("-a","--action",required=True,type=str,choices=["import-smiles","pick","import-dock","train","sample-pred","pred","import-pred"],help="Action to perform")
    parser.add_argument("-c","--cpu",required=False,type=int,help="How many CPUs to use")
    parser.add_argument("-s","--smiles",required=False,type=str,help="SMILES filename when importing")
    parser.add_argument("-n","--name",required=False,type=str,help="Name for the virtual screening")
    parser.add_argument("-f","--fraction",required=False,type=float,help="Fraction of database to dock in each iteration")
    parser.add_argument("-d","--dockingtemplate",required=False,type=str,help="Filename of docking template")
    parser.add_argument("-x","--predchunksize",required=False,type=int,default=123456,help="Size of chemprop prediction chunk")
    return parser.parse_args()

def files_exist(args):
    """
    Check if input files exist

    :param args: parsed arguments
    :return: True if files exist, false otherwise
    """
    if args.database is not None and not os.path.exists(args.database) and args.action!="import-smiles":
        print("HASTEN database file missing!")
        return False
    if args.action=="import-smiles" and (args.smiles==None or not os.path.exists(args.smiles)):
        print("SMILES file missing while trying import to database")
        return False
    if args.action=="train" and (args.iteration==None or args.name==None or not os.path.exists("hasten_"+args.name+".db")):
        print("Error in command line while trying to produce training input. Check that you have all parameters and that name is correct.")
        return False
    return True

def pick_compounds_for_docking(db,iteration,name,fraction,cpu,dockingtemplate):
    """
    Pick set of compounds from database for docking and save them into
    output files.

    :param db: The filename of SQlite3 database
    :param iteration: iteration integer
    :param name: name of the virtual screen
    :param fraction: fraction of molecules to dock
    :param cpu: number of CPUs to use
    :param dockingtemplate: template for docking
    """
    if name == None or fraction == None or cpu == None:
        print("Error while picking compounds for docking: check your commmand line.")
        sys.exit(1)
    conn=sqlite3.connect(db)
    c = conn.cursor()
    # counting is very slow and we are not deleting so use a hack here
    # https://stackoverflow.com/questions/8988915/sqlite-count-slow-on-big-tables
    #number_of_mols=c.execute("SELECT COUNT() FROM data").fetchone()[0]
    number_of_mols=c.execute("SELECT MAX(_ROWID_) FROM data LIMIT 1").fetchone()[0]
    number_to_dock=int(round(fraction*number_of_mols))
    print("Number of molecules in the database",number_of_mols)
    print("Picking",fraction*100,"% for docking (",number_to_dock,")...")
    number_per_cpu = int(number_to_dock/cpu) + (number_to_dock % cpu >0)
    print("Molecules per CPU:",number_per_cpu)

    # first time pick just random set (slow way but as this is done once in
    # every screen only it does not matter)
    if iteration==1:
        to_dock=c.execute("SELECT smiles,hastenid FROM data ORDER BY random() LIMIT ?",[number_to_dock,]).fetchall()
    else:
        if not os.path.exists("hasten_"+str(name)+".db") or not os.path.exists("predictions.db"):
            print("Error: docking results and/or predictions database(s) missing.")
            sys.exit(1)
        c.execute('ATTACH DATABASE "hasten_'+str(name)+'.db" AS d')
        c.execute('ATTACH DATABASE "predictions.db" AS p')
        to_dock=c.execute("SELECT smiles,data.hastenid FROM data,pred_data WHERE data.hastenid=pred_data.hastenid AND data.hastenid NOT IN (SELECT hastenid FROM docking_data) ORDER BY pred_data.pred_score LIMIT ?",[number_to_dock,]).fetchall()

    w = open(name+"_iter"+str(iteration)+"_cpu1.smi","wt")
    # write smilesids
    mols_per_file = 0
    cpu_counter = 1
    while (len(to_dock) > 0):
        mols_per_file += 1
        smiles,hastenid=to_dock.pop()
        w.write(str(smiles)+" "+str(hastenid)+"\n")
        if mols_per_file >= number_per_cpu:
            cpu_counter +=1
            mols_per_file = 0
            w.close()
            if len(to_dock)>0:
                w = open(name+"_iter"+str(iteration)+"_cpu"+str(cpu_counter)+".smi","wt")        
    w.close()

    # create inputs for conformation generation
    for i in range(1,cpu+1):
        w = open(name+"_iter"+str(iteration)+"_cpu"+str(i)+".inp","wt")
        w.write("[SET:ORIGINAL_LIGANDS]\n")
        w.write("    VARCLASS   Structures\n")
        w.write("    FILES   "+name+"_iter"+str(iteration)+"_cpu"+str(i)+".smi"+",\n")
        w.write("\n")
        w.write("[STAGE:LIGPREP]\n")
        w.write("    STAGECLASS   ligprep.LigPrepStage\n")
        w.write("    INPUTS   ORIGINAL_LIGANDS,\n")
        w.write("    OUTPUTS   LIGPREP_OUT,\n")
        w.write("    RECOMBINE   YES\n")
        w.write("    RETITLE   YES\n")
        w.write("    MIXLIGS   YES\n")
        w.write("    SKIP_BAD_LIGANDS   YES\n")
        w.write("    UNIQUEFIELD   s_m_title\n")
        w.write("    OUTCOMPOUNDFIELD   s_m_title\n")
        w.write("    USE_EPIK   YES\n")
        w.write("    METAL_BINDING   NO\n")
        w.write("    PH   7.0\n")
        w.write("    PHT   2.0\n")
        w.write("    NRINGCONFS   1\n")
        w.write("    COMBINEOUTS   NO\n")
        w.write("    STEREO_SOURCE   parities\n")
        w.write("    NUM_STEREOISOMERS   32\n")
        w.write("    REGULARIZE   NO\n")
        w.write("\n")
        w.write("[STAGE:POSTLIGPREP]\n")
        w.write("    STAGECLASS   ligprep.PostLigPrepStage\n")
        w.write("    INPUTS   LIGPREP_OUT,\n")
        w.write("    OUTPUTS   POSTLIGPREP_OUT,\n")
        w.write("    UNIQUEFIELD   s_m_title\n")
        w.write("    OUTVARIANTFIELD   s_phase_variant\n")
        w.write("    PRESERVE_NJOBS   YES\n")
        w.write("    LIMIT_STEREOISOMERS   YES\n")
        w.write("    MAXSTEREO   4\n")
        w.write("    REMOVE_PENALIZED_STATES   YES\n")
        w.write("\n")
        w.write("[STAGE:MANAGE]\n")
        w.write("    STAGECLASS   phase.DBManageStage\n")
        w.write("    INPUTS   POSTLIGPREP_OUT,\n")
        w.write("    OUTPUTS   DATABASE,\n")
        w.write("    DATABASE  "+name+"_iter"+str(iteration)+"_cpu"+str(i)+".phdb"+"\n")
        w.write("    NEW   YES\n")
        w.write("    MULTIPLE_CONFS   NO\n")
        w.write("    CONSIDER_STEREO   NO\n")
        w.write("    GENERATE_PROPS   NO\n")
        w.write("    CREATE_SUBSET   NO\n")
        w.write("    SKIP_DUPLICATES   NO\n")
        w.write("\n")
        w.write("[STAGE:CONFSITES]\n")
        w.write("    STAGECLASS   phase.DBConfSitesStage\n")
        w.write("    INPUTS   DATABASE,\n")
        w.write("    CONFS   auto\n")
        w.write("    MAX_CONFS   1\n")
        w.write("    GENERATE_PROPS   YES\n")
        w.write("\n")
        w.write("[USEROUTS]\n")
        w.write("    USEROUTS   DATABASE,\n")
        w.close()

    template_lines = []
    for line in open(dockingtemplate,"rt"):
        l = line.strip().split()
        if len(l)>0 and l[0] != "LIGANDFILE":
            template_lines.append(line)
    # create inputs for docking
    for i in range(1,cpu+1):
        w = open("glide_"+name+"_iter"+str(iteration)+"_cpu"+str(i)+".in","wt")
        w.write("LIGANDFILE   "+name+"_iter"+str(iteration)+"_cpu"+str(i)+"_1.maegz\n")
        for template_line in template_lines:
            w.write(template_line)
        w.close()

def import_smiles(database,smilesfile):
    print("Importing",smilesfile,end="")
    if not os.path.exists(database):
        print(" to a new database ",end="")
    else:
        print(" to a existing database ",end="")
    print(database)
    conn=sqlite3.connect(args.database)
    c=conn.cursor()
    c.execute("CREATE TABLE IF NOT EXISTS data (hastenid INTEGER PRIMARY KEY,smiles TEXT,smilesid TEXT)")
    chunksize=1000000
    to_db = []
    imported=0
    with open(args.smiles) as smilesfile:
        for row in csv.reader(smilesfile,delimiter=" "):
            to_db.append((row[0],row[1]))
            if len(to_db)>=chunksize:
                c.executemany("INSERT INTO data(smiles,smilesid) VALUES (?,?)",to_db)
                imported+=len(to_db)
                to_db = []
                print(imported,"molecules imported...")
        conn.commit()
    if len(to_db)>0:
        c.executemany("INSERT INTO data(smiles,smilesid) VALUES (?,?)",to_db)
        conn.commit()
        imported+=len(to_db)
    conn.close()
    print("Done, imported",imported,"molecules")

def import_dock(database,iteration,name):
    conn = sqlite3.connect("hasten_"+name+".db")
    c = conn.cursor()
    c.execute("CREATE TABLE IF NOT EXISTS docking_data (hastenid INTEGER PRIMARY KEY,dock_score NUMERIC,dock_iteration INTEGER)")
    conn.commit()

    print("Importing results to hasten_"+name+".db...")
    docking_scores = {}
    for filename in glob.glob("glide_"+name+"_iter"+str(iteration)+"_cpu*.csv"):
        with open(filename) as smilesfile:
            for row in csv.DictReader(smilesfile,delimiter=","):
                hastenid = int(row["title"])
                try:
                    ds = float(row["r_i_docking_score"])
                    if ds > 10.0:
                        ds = 10.0
                except:
                    ds = 10.0
                if hastenid not in docking_scores or docking_scores[hastenid]>ds:
                    docking_scores[hastenid] = ds
    to_db = []
    for hastenid in docking_scores:
        to_db.append((hastenid,docking_scores[hastenid],iteration))
    c.executemany("INSERT OR REPLACE INTO docking_data(hastenid,dock_score,dock_iteration) VALUES (?,?,?)",to_db)
    conn.commit()

def import_pred(name):
    conn = sqlite3.connect("predictions.db")
    c = conn.cursor()
    c.execute("CREATE TABLE IF NOT EXISTS pred_data (hastenid INTEGER PRIMARY KEY,pred_score NUMERIC)")
    conn.commit()

    count = 0
    print("Importing predictions to hasten_"+name+".db...")
    for predpath in glob.glob("PRED*"):
        count += 1
        print(count,predpath)
        for filename in glob.glob(predpath+"/output_hasten_pred_tmp*.csv"):
            predicted_scores = []
            with open(filename) as smilesfile:
                for row in csv.DictReader(smilesfile,delimiter=","):
                    predicted_scores.append((int(row["hastenid"]),float(row["docking_score"])))
            c.executemany("INSERT OR REPLACE INTO pred_data(hastenid,pred_score) VALUES (?,?)",predicted_scores)
            conn.commit()

def train(database,iteration,name):
    s_conn = sqlite3.connect(database)
    sc = s_conn.cursor()
    sc.execute('ATTACH DATABASE "hasten_'+str(name)+'.db" AS d')
    s_conn.commit()
    smilesdata=sc.execute("SELECT smiles,data.hastenid,docking_data.dock_score FROM data,docking_data WHERE dock_score IS NOT NULL AND data.hastenid=docking_data.hastenid").fetchall()
    w = open("train_"+name+"_iter"+str(iteration)+".csv","wt")
    w.write("smiles,hastenid,docking_score\n")
    for smiles,hastenid,docking_score in smilesdata:
        w.write(smiles+","+str(hastenid)+","+str(docking_score)+"\n")
    w.close()

def pred(database,iteration,cpu,name,predchunksize):
    if not os.path.exists("cutoff_"+name+"_iter"+str(iteration)+".txt"):
        print("Cutoff file missing, run sample-pred before this step.")
        sys.exit(1)
    cutoff = open("cutoff_"+name+"_iter"+str(iteration)+".txt").readlines()[0].strip()
    conn = sqlite3.connect(database)
    c = conn.cursor()
    number_of_mols=c.execute("SELECT MAX(_ROWID_) FROM data LIMIT 1").fetchone()[0]
    print("Number of molecules in the database",number_of_mols)
    number_per_cpu = int(number_of_mols/cpu) + (number_of_mols % cpu >0)
    print("Molecules per CPU:",number_per_cpu)

    mols_left = number_of_mols
    chunk_size = predchunksize
    cur_index = 1

    jobs = []
    cpu_counter = 1
    while cur_index < mols_left:
        jobs.append((cpu_counter,'sqlite3 -separator "," '+database+' "SELECT smiles,hastenid FROM data WHERE hastenid>='+str(cur_index)+' AND hastenid<'+str(cur_index+chunk_size)+'"'))
        cur_index += chunk_size
        cpu_counter += 1
        if cpu_counter > cpu:
            cpu_counter = 1

    cur_cpu = 1
    w = open("pred_"+name+"_iter"+str(iteration)+"_cpu"+str(cur_cpu)+".sh","wt")
    if not os.path.exists("PRED"+str(cur_cpu)):
        os.mkdir("PRED"+str(cur_cpu))
    tmp_counter = 0
        
    for scpu,job in sorted(jobs):
        tmp_counter +=1
        if scpu>cur_cpu:
            w.close()
            cur_cpu+=1
            w = open("pred_"+name+"_iter"+str(iteration)+"_cpu"+str(cur_cpu)+".sh","wt")
            if not os.path.exists("PRED"+str(cur_cpu)):
                os.mkdir("PRED"+str(cur_cpu))
        tmpfile="hasten_pred_tmp"+str(tmp_counter)+".csv"
        w.write("echo smiles,hastenid >PRED"+str(scpu)+"/"+tmpfile+"\n")
        w.write(job+">>PRED"+str(scpu)+"/"+tmpfile+"\n")
        w.write("OMP_NUM_THREADS=1 nice chemprop_predict --num_workers 0 --no_cache_mol --no_cuda --test_path PRED"+str(scpu)+"/"+tmpfile+" --checkpoint_dir "+name+"_iter"+str(iteration)+" --preds_path PRED"+str(scpu)+"/output_"+tmpfile+"\n")
        w.write("awk -F, '{ if ($3<="+str(cutoff)+") { print } }' PRED"+str(scpu)+"/output_"+tmpfile+" >PRED"+str(scpu)+"/tmp_"+tmpfile+"\n")
        w.write("echo smiles,hastenid,docking_score >PRED"+str(scpu)+"/output_"+tmpfile+"\n")
        w.write("cat PRED"+str(scpu)+"/tmp_"+tmpfile+" >>PRED"+str(scpu)+"/output_"+tmpfile+"\n")
        w.write("rm PRED"+str(scpu)+"/"+tmpfile+" PRED"+str(scpu)+"/tmp_"+tmpfile+"\n")
    w.close()

def sample_pred(database,iteration,name,cpu):
    if os.path.exists("output_samplepred_"+name+"_iter"+str(iteration)+"_cpu1.csv"):
        print("Processing sampled predictions...")
        pred_scores = []
        skipped = 0
        for filename in glob.glob("output_samplepred_"+name+"_iter"+str(iteration)+"_cpu*.csv"):
            with open(filename) as outputfile:
                csvreader = csv.DictReader(outputfile,delimiter=",")
                for row in csvreader:
                    try:
                        pred_scores.append(float(row["docking_score"]))
                    except:
                        print("NOTE: skipped",row)
                        skipped += 1
        x = statistics.mean(pred_scores)
        sd = statistics.stdev(pred_scores)
        cutoff = round(x-0.5*sd,1)
        print("Scores used:",len(pred_scores))
        print("Skipped SMILES:",skipped)
        print("Mean of predicted docking scores:",x)
        print("Sample standard deviation of data:",sd)
        print("\nPrediction cutoff:",cutoff)
        w = open("cutoff_"+name+"_iter"+str(iteration)+".txt","wt")
        w.write(str(cutoff)+"\n")
        w.close()
    else:
        conn=sqlite3.connect(database)
        c = conn.cursor()

        number_of_mols=c.execute("SELECT MAX(_ROWID_) FROM data LIMIT 1").fetchone()[0]
        number_to_pred=int(round(0.001*number_of_mols))
        print("Picking",number_to_pred,"for prediction cutoff sampling...")
        number_per_cpu = int(number_to_pred/cpu) + (number_to_pred % cpu >0)

        to_pred=c.execute("SELECT smiles,hastenid FROM data ORDER BY random() LIMIT ?",[number_to_pred,]).fetchall()

        w = open("samplepred_"+name+"_iter"+str(iteration)+"_cpu1.csv","wt")        
        w.write("smiles,hastenid\n")
        # write smilesids
        mols_per_file = 0
        cpu_counter = 1
        while (len(to_pred) > 0):
            mols_per_file += 1
            smiles,hastenid=to_pred.pop()
            w.write(str(smiles)+","+str(hastenid)+"\n")
            if mols_per_file >= number_per_cpu:
                cpu_counter +=1
                mols_per_file = 0
                w.close()
                if len(to_pred)>0:
                    w = open("samplepred_"+name+"_iter"+str(iteration)+"_cpu"+str(cpu_counter)+".csv","wt")        
                    w.write("smiles,hastenid\n")
        w.close()
 
def run_hasten(args):
    """
    Starts the whole process

    :param args: Parsed arguments
    """
    if args.iteration is not None:
        iteration = args.iteration
    else:
        iteration = 1

    print("Iteration",iteration)
    if args.action == "import-smiles":
        import_smiles(args.database,args.smiles)
    elif args.action == "pick":
        pick_compounds_for_docking(args.database,args.iteration,args.name,args.fraction,args.cpu,args.dockingtemplate)
    elif args.action == "import-dock":
        import_dock(args.database,args.iteration,args.name)
    elif args.action == "train":
        train(args.database,args.iteration,args.name)
    elif args.action == "pred":
        pred(args.database,args.iteration,args.cpu,args.name,args.predchunksize)
    elif args.action == "sample-pred":
        sample_pred(args.database,args.iteration,args.name,args.cpu)
    elif args.action == "import-pred":
        import_pred(args.name)
    else:
        print("ERROR IN CODE: UNDEFINED ACTION ALLOWED.")
        sys.exit(1)

    print("\nHASTEN finished.")

if __name__ == "__main__":
    print("")
    print("    )        (                 )  ")
    print(" ( /(  (     )\ )  *   )    ( /(  ")
    print(" )\()) )\   (()/(` )  /((   )\()) ")
    print("((_)((((_)(  /(_))( )(_))\ ((_)\  ")
    print(" _((_)\ _ )\(_)) (_(_()|(_) _((_) ")
    print("| || (_)_\(_) __||_   _| __| \| | ")
    print("| __ |/ _ \ \__ \  | | | _|| .` | ")
    print("|_||_/_/ \_\|___/  |_| |___|_|\_| ")
    print("")
    print("HASTEN (macHine leArning booSTEd dockiNg) version 1.0")
    print("https://github.com/TuomoKalliokoski/HASTEN")
    print("")
    print("Reference:")
    print("Kalliokoski T. Molecular Informatics 2021, doi:10.1002/minf.202100089")
    print("")
    args = parse_cmd_line()
    if (not files_exist(args)): sys.exit(1)
    run_hasten(args)

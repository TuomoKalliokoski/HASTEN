# 
# Copyright (c) 2021-2023 Orion Corporation
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

v1.1.2:
    - support for importing gzip compressed SMILES added
    - fixed couple of variable names in import_smiles (reported by lewisri2)

v1.1:
    - failed dockings are excluded from the training [study by Sivula et al]. Switch "-t" added
    - better on-line help implemented
    - predictions are now imported into seperate databases each run
    - progress counter added to the prediction importing
    - added status functionality that allows checking the number of hits quickly
    - add exporting of results both in SMILES and MAEGZ formats
    - import SMILES now checks that the input has at least two space separated fields
    - SMILES are removed from prediction output (saves disk space)
    - gawk forced

"""

hastentutorial="""

You need three input files before you can start HASTEN:
    1. database to dock in SMILES (space delimited) (can be gzip compressed)
    2. Glide docking template (.in file)
    3. Glide docking grid that is defined in #2

Start the HASTEN procedure by importing the SMILES to HASTEN DB:
python3 hasten.py -m mols.db -a import-smiles -s mols.smi

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
import gzip

def parse_cmd_line():
    """
    Parse command line using ArgumentParser

    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description=hastentutorial)
    parser.add_argument("-m","--database",required=True,type=str,help="HASTEN database")
    parser.add_argument("-i","--iteration",required=False,type=int,help="Iteration number")

    parser.add_argument("-a","--action",required=True,type=str,choices=["import-smiles","pick","import-dock","train","sample-pred","pred","import-pred","status","export-smiles","export-poses"],help="Action to perform")
    parser.add_argument("-c","--cpu",default=1,required=False,type=int,help="How many CPUs to use")
    parser.add_argument("-s","--smiles",required=False,type=str,help="SMILES filename when importing/exporting (.gz supported in import)")
    parser.add_argument("-n","--name",required=False,type=str,help="Name for the virtual screening")
    parser.add_argument("-f","--fraction",required=False,type=float,help="Fraction of database to dock in each iteration")
    parser.add_argument("-d","--dockingtemplate",required=False,type=str,help="Filename of docking template")
    parser.add_argument("-x","--predchunksize",required=False,type=int,default=123456,help="Size of chemprop prediction chunk")
    parser.add_argument("-t","--failed",required=False,type=float,default=None,help="(Optional) By default, failed dockings are NOT included in training. You can set include them with desired value using this switch")
    parser.add_argument("-u","--cutoff",required=False,type=float,default=None,help="(Optional) Set exporting cutoff (otherwise just virtual hits are exported)")
    args = parser.parse_args()
    # make sure that we keep the 1909 value reserved for failed dockings
    if args.failed and args.failed>=12345.0:
        parse.error("Too large failed score (must be below 12345.0)")
    return args

def assist_user(short_help,help_string):
    print("\nNext step ("+short_help+"):")
    print(help_string)

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
        pred_dbname = "predictions_"+str(name)+"_"+str(iteration)+".db"
        if not os.path.exists("hasten_"+str(name)+".db") or not os.path.exists(pred_dbname):
            print("Error: docking (hasten_"+str(name)+".db) and/or predictions database(s) ("+pred_dbname+") missing.")
            sys.exit(1)
        c.execute('ATTACH DATABASE "hasten_'+str(name)+'.db" AS d')
        c.execute('ATTACH DATABASE "'+pred_dbname+'" AS p')
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
    assist_user("Phase LigPrep","find "+name+"_iter"+str(iteration)+"_*.inp -exec echo nice $SCHRODINGER/pipeline -prog phase_db {} -OVERWRITE -WAIT -HOST localhost:1 -NJOBS 1 \; | parallel --bar -j "+str(cpu))
    assist_user("convert to maegz","find `pwd` -maxdepth 1 -name \""+name+"_iter"+str(iteration)+"_*.phdb\" -exec echo $SCHRODINGER/phase_database {} export -omae {} \; | sed -e 's:.phdb$: -get 1 -limit 99999999 -WAIT:' | parallel --bar -j "+str(cpu))
    assist_user("Glide","find glide_"+name+"_iter"+str(iteration)+"_*.in -exec echo nice $SCHRODINGER/glide {} -NICE -OVERWRITE -WAIT -NJOBS 1 -HOST localhost:1 \; | parallel --bar -j "+str(cpu))
    assist_user("Import docking results","python3 "+sys.argv[0]+" -m "+db+" -a import-dock -i "+str(iteration)+" -n "+name)

def import_smiles(database,smilesfile):
    """
    Import SMILES to HASTEN database

    :param database: .db filename where to compounds are to be imported
    :param smilesfile: .smi filename to be imported
    """
    print("Importing",smilesfile,end="")
    if not os.path.exists(database):
        print(" to a new database ",end="")
    else:
        print(" to an existing database ",end="")
    print(database)
    conn=sqlite3.connect(database)
    c=conn.cursor()
    c.execute("CREATE TABLE IF NOT EXISTS data (hastenid INTEGER PRIMARY KEY,smiles TEXT,smilesid TEXT)")
    chunksize=1000000
    to_db = []
    imported=0
    opener = open
    if ".gz" in smilesfile:
        opener = gzip.open
    with opener(args.smiles,"rt") as smilesfile:
        for row in csv.reader(smilesfile,delimiter=" "):
            if len(row)<2:
                print("Error: the input SMILES file should be space-delimited file, import stopped.")
                sys.exit(1)
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
    assist_user("pick 1st iteration compounds to be docked","python3 "+sys.argv[0]+" -m "+database+" -a pick -i 1 -n <jobname> -f <fraction_to_dock> -c <cpus_to_use> -d <template.in>")
    
def import_dock(database,iteration,name,failed):
    """
    Import docking results

    :param database: .db filename where to compounds are to be imported
    :param iteration: iteration integer
    :param name: name of the virtual screen
    :param failed: how to handle failed dockings (None=default whic is to ignore and store as 1909.0 for failed)
    """
    conn = sqlite3.connect("hasten_"+name+".db")
    c = conn.cursor()
    c.execute("CREATE TABLE IF NOT EXISTS docking_data (hastenid INTEGER PRIMARY KEY,dock_score NUMERIC,dock_iteration INTEGER)")
    conn.commit()

    docking_scores = {}
    filenames = glob.glob("glide_"+name+"_iter"+str(iteration)+"_cpu*.csv")
    file_counter = 0
    for filename in filenames:
        file_counter +=1
        print("Importing results to hasten_"+name+".db,",file_counter,"of",len(filenames),":",filename)
        with open(filename) as smilesfile:
            for row in csv.DictReader(smilesfile,delimiter=","):
                hastenid = int(row["title"])
                try:
                    ds = float(row["r_i_docking_score"])
                    # clamp the docking values
                    if failed!=None and ds>failed:
                        ds = failed
                except:
                    if failed!=None:
                        ds = failed
                    else:
                        # use this constant to mark compound as failed
                        ds = 12345.0
                if hastenid not in docking_scores or docking_scores[hastenid]>ds:
                    docking_scores[hastenid] = ds
    to_db = []
    for hastenid in docking_scores:
        to_db.append((hastenid,docking_scores[hastenid],iteration))
    c.executemany("INSERT OR REPLACE INTO docking_data(hastenid,dock_score,dock_iteration) VALUES (?,?,?)",to_db)
    conn.commit()
    assist_user("quick status","python3 "+sys.argv[0]+" -m "+database+" -a status -n "+name)
    assist_user("Go to iteration "+str(iteration+1),"python3 "+sys.argv[0]+" -a train -m "+database+" -n "+name+" -i "+str(iteration+1))

def import_pred(database,name,iteration):
    """
    Import ML prediction results

    :param database: the name of the main HASTEN database for the assist_user
    :param name: name of the virtual screen
    :param iteration: iteration
    """
    pred_dbname="predictions_"+name+"_"+str(iteration)+".db"
    if os.path.exists(pred_dbname):
        os.remove(pred_dbname)
    conn = sqlite3.connect(pred_dbname)
    c = conn.cursor()
    c.execute("CREATE TABLE pred_data (hastenid INTEGER PRIMARY KEY,pred_score NUMERIC)")
    conn.commit()

    count = 0
    print("Importing predictions to "+pred_dbname+"...")
    res_dirs = glob.glob("PRED*")
    for predpath in res_dirs:
        count += 1
        predfiles=glob.glob(predpath+"/output_hasten_pred_tmp*.csv")
        icount = 0 
        for filename in predfiles:
            icount += 1
            print("Processing:",predpath,"[",count,"of",len(res_dirs),"][",icount,"of",len(predfiles),"]")
            predicted_scores = []
            with open(filename) as smilesfile:
                for row in csv.DictReader(smilesfile,delimiter=","):
                    predicted_scores.append((int(row["hastenid"]),float(row["docking_score"])))
            c.executemany("INSERT INTO pred_data(hastenid,pred_score) VALUES (?,?)",predicted_scores)
            conn.commit()
    assist_user("Pick compounds for docking","python3 "+sys.argv[0]+" -m "+database+" -a pick -i "+str(iteration)+" -n "+name+" -f <fraction_to_dock> -c <cpu_count> -d <template.in>")

def train(database,iteration,name,failed):
    """
    Train ML model

    :param database: The filename of SQlite3 database
    :param iteration: iteration integer
    :param name: name of the virtual screen
    :param failed: how to deal with failed scores (None = don't use in the training)
    """
    print("Preparing input file for ML training...")
    s_conn = sqlite3.connect(database)
    sc = s_conn.cursor()
    sc.execute('ATTACH DATABASE "hasten_'+str(name)+'.db" AS d')
    s_conn.commit()
    if failed==None:
        print("Failed compounds not included")
    else:
        print("Failed compounds included at cutoff",failed)
    if failed==None:
        smilesdata=sc.execute("SELECT smiles,data.hastenid,docking_data.dock_score FROM data,docking_data WHERE dock_score IS NOT NULL AND dock_score < 1909 AND data.hastenid=docking_data.hastenid").fetchall()
    else:
        smilesdata=sc.execute("SELECT smiles,data.hastenid,docking_data.dock_score FROM data,docking_data WHERE dock_score IS NOT NULL AND data.hastenid=docking_data.hastenid").fetchall()
    w = open("train_"+name+"_iter"+str(iteration)+".csv","wt")
    w.write("smiles,hastenid,docking_score\n")
    for smiles,hastenid,docking_score in smilesdata:
        w.write(smiles+","+str(hastenid)+","+str(docking_score)+"\n")
    w.close()
    assist_user("start ChemProp training","chemprop_train --dataset_type regression --target_columns docking_score --data_path train_"+name+"_iter"+str(iteration)+".csv --save_dir "+name+"_iter"+str(iteration)+" --batch_size 250 --no_cache_mol")
    assist_user("pick sample predictions","python3 "+sys.argv[0]+" -m "+database+" -a sample-pred -i "+str(iteration)+" -n "+name+" -c 1")

def pred(database,iteration,cpu,name,predchunksize):
    """
    Predict scores using ML model

    :param database: The filename of SQlite3 database
    :param iteration: iteration integer
    :param cpu: Number of cores
    :param name: name of the virtual screen
    :param predchunksize: how many compounds per prediction chunk
    """
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
        w.write("OMP_NUM_THREADS=1 nice chemprop_predict --num_workers 0 --no_cache_mol --no_cuda --test_path PRED"+str(scpu)+"/"+tmpfile+" --checkpoint_dir "+name+"_iter"+str(iteration)+" --preds_path PRED"+str(scpu)+"/tmp_"+tmpfile+"\n")
        w.write("gawk -F, 'NR==1 { print $2\",\"$3 } NR>1 { if ($3<="+str(cutoff)+") { print $2\",\"$3 } }' PRED"+str(scpu)+"/tmp_"+tmpfile+" >PRED"+str(scpu)+"/output_"+tmpfile+"\n")
        w.write("rm PRED"+str(scpu)+"/"+tmpfile+" PRED"+str(scpu)+"/tmp_"+tmpfile+"\n")
    w.close()
    assist_user("Make sure you have gawk","sudo apt install gawk")
    assist_user("If you want to use GPU instead of CPU","sed -i -e 's:OMP_NUM_THREADS=1::' -e 's:--num_workers 0::' -e 's:--no_cache_mol --no_cuda::' pred_"+name+"_iter"+str(iteration)+"_cpu*.sh")
    assist_user("run ChemProp predict","ls -1 pred_"+name+"_iter"+str(iteration)+"_cpu*.sh | sed -e 's:^:sh :' | parallel -j "+str(cpu)+" --bar")
    assist_user("import predictions","python3 "+sys.argv[0]+" -a import-pred -i "+str(iteration)+" -n "+name)

def sample_pred(database,iteration,name,cpu):
    """
    Do sample prediction in order to set the cutoff for throwing predictions away

    :param database: The filename of SQlite3 database
    :param iteration: iteration integer
    :param name: name of the virtual screen
    :param cpu: number of CPUs to use
    """
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
        assist_user("prepare to predict","python3 "+sys.argv[0]+" -m "+database+" -a pred -i "+str(iteration)+" -n "+name+" -c <cpu_count> -x 10000")
    else:
        conn=sqlite3.connect(database)
        c = conn.cursor()

        number_of_mols=c.execute("SELECT MAX(_ROWID_) FROM data LIMIT 1").fetchone()[0]
        number_to_pred=int(round(0.001*number_of_mols))
        print("Picking",number_to_pred,"for prediction cutoff sampling...")
        assist_user("run sample predictions","ls -1 samplepred_"+name+"_iter"+str(iteration)+"_cpu*.csv | awk \'{ print \"chemprop_predict --checkpoint_dir "+name+"_iter"+str(iteration)+" --test_path\",$1,\"--preds_path output_\"$1 }\' | parallel --bar -j 1")
        assist_user("set prediction cutoff","python3 "+sys.argv[0]+" -m "+database+" -a sample-pred -i "+str(iteration)+" -n "+name+" -c 1")
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
 
def show_status(database,name):
    docking_db = "hasten_"+str(name)+".db"
    if not os.path.exists(docking_db):
        print("No docking results found for this job ("+docking_db+")")
    print("Virtual screening jobname:",name)
    print("Docking scores database:",docking_db)
    conn = sqlite3.connect("hasten_"+name+".db")
    c = conn.cursor()
    max_iteration = c.execute("SELECT MAX(dock_iteration) FROM docking_data LIMIT 1").fetchone()[0]
    best_docking_score = c.execute("SELECT MIN(dock_score) FROM docking_data LIMIT 1").fetchone()[0]
    virtual_hit = 0.85 * best_docking_score
    print("\nDocking iterations done:",max_iteration)
    print("Best docking score:",best_docking_score)
    print("Virtual hit cutoff",virtual_hit)
    print("\nIter\tVirtual hits")
    all_hits = 0
    for i in range(max_iteration):
        virtual_hits = c.execute("SELECT COUNT(*) FROM docking_data WHERE dock_iteration LIKE "+str(i+1)+" AND dock_score <= "+str(virtual_hit)).fetchone()[0]
        all_hits += virtual_hits
        print(str(i+1)+"\t"+str(virtual_hits))
    print("\nNumber of virtual hits found:",all_hits)
    assist_user("Export virtual hits to SMILES","python3 "+sys.argv[0]+" -a export-smiles -m "+database+" -n "+name+" -s virtualhits_"+name+".smi -u "+str(virtual_hit))
    c.close()
    conn.close()
    
def export_smiles(database,name,smiles,cutoff):
    s_conn = sqlite3.connect(database)
    sc = s_conn.cursor()
    sc.execute('ATTACH DATABASE "hasten_'+str(name)+'.db" AS d')
    s_conn.commit()
    if cutoff==None:
        best_docking_score = sc.execute("SELECT MIN(dock_score) FROM docking_data LIMIT 1").fetchone()[0]
        virtual_hit = 0.85 * best_docking_score
    else:
        virtual_hit = cutoff
    print("Docking score cut-off for exporting:",virtual_hit)
    print("Output file:",smiles)
    smilesdata=sc.execute("SELECT smiles,data.smilesid,docking_data.dock_score,docking_data.dock_iteration FROM data,docking_data WHERE dock_score <= "+str(virtual_hit)+" AND data.hastenid=docking_data.hastenid").fetchall()
    w = open(smiles,"wt")
    for smiles,smilesid,docking_score,dock_iteration in smilesdata:
        w.write(smiles+" "+str(smilesid)+" "+str(docking_score)+" "+str(dock_iteration)+"\n")
    w.close()
    assist_user("export docking poses","$SCHRODINGER/run "+sys.argv[0]+" -m "+database+" -a export-poses -n "+name+" -u "+str(virtual_hit))

def export_poses(database,name,cutoff):
    try:
        from schrodinger import structure
    except:
        print("Schrödinger Suite Python environment not detected, did you start HASTEN with $SCHRODINGER/run?")
        sys.exit(1)
    s_conn = sqlite3.connect(database)
    sc = s_conn.cursor()
    sc.execute('ATTACH DATABASE "hasten_'+str(name)+'.db" AS d')
    s_conn.commit()
    if cutoff==None:
        best_docking_score = sc.execute("SELECT MIN(dock_score) FROM docking_data LIMIT 1").fetchone()[0]
        virtual_hit = 0.85 * best_docking_score
    else:
        virtual_hit = cutoff
    print("Docking score cut-off for exporting:",virtual_hit)
    output_filename = "virtualhits_"+name+".maegz"
    print("Output file:",output_filename)
    writer = structure.StructureWriter(output_filename)
    resfiles = glob.glob("glide_"+name+"_iter*.maegz")

    hit_ids = set()
    poses = 0
    count = 0
    for resfile in resfiles:
        count += 1
        print("Processing",count,"of",len(resfiles),":",resfile)
        reader = structure.StructureReader(resfile)
        first = True
        for st in reader:
            # skip protein as the first structure if pv file was produced
            if "pv" in resfile and first:
                first = False
                continue
            code_id = int(st.property["s_m_title"])
            docking_score = float(st.property["r_i_docking_score"])
            if docking_score <= virtual_hit:
                st.property["s_m_title"] = sc.execute("SELECT smilesid FROM data WHERE hastenid="+str(code_id)).fetchone()[0]
                writer.append(st)
                poses += 1
                if st.property["s_m_title"] not in hit_ids:
                    hit_ids.add(st.property["s_m_title"])
            elif "pv" in resfile:
                # pv files are sorted so we can get out, otherwise we need to go through the file (raw files without sorting)
                break
            
    writer.close()
    s_conn.close()
    print("Virtual hits written:",len(hit_ids),"("+str(poses)+" poses)")

def run_hasten(args):
    """
    Starts the whole process

    :param args: Parsed arguments
    """
    if args.iteration is not None:
        iteration = args.iteration
    else:
        iteration = 1

    if args.action != "import-smiles" and args.action!="status" and args.action!="export-smiles" and args.action!="export-poses":
        print("Iteration",iteration)

    if args.action == "import-smiles":
        import_smiles(args.database,args.smiles)
    elif args.action == "pick":
        pick_compounds_for_docking(args.database,args.iteration,args.name,args.fraction,args.cpu,args.dockingtemplate)
    elif args.action == "import-dock":
        import_dock(args.database,args.iteration,args.name,args.failed)
    elif args.action == "train":
        train(args.database,args.iteration,args.name,args.failed)
    elif args.action == "pred":
        pred(args.database,args.iteration,args.cpu,args.name,args.predchunksize)
    elif args.action == "sample-pred":
        sample_pred(args.database,args.iteration,args.name,args.cpu)
    elif args.action == "import-pred":
        import_pred(args.database,args.name,args.iteration)
    elif args.action == "status":
        show_status(args.database,args.name)
    elif args.action == "export-smiles":
        export_smiles(args.database,args.name,args.smiles,args.cutoff)
    elif args.action == "export-poses":
        export_poses(args.database,args.name,args.cutoff)
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
    print("HASTEN (macHine leArning booSTEd dockiNg) version 1.1.2")
    print("https://github.com/TuomoKalliokoski/HASTEN")
    print("")
    print("References:")
    print("Kalliokoski T. Mol. Inform. 2021, doi:10.1002/minf.202100089")
    print("Sivula T et al.  J. Chem. Inf. Model. 2023, doi:10.1021/acs.jcim.3c01239")
    print("")
    args = parse_cmd_line()
    if (not files_exist(args)): sys.exit(1)
    run_hasten(args)

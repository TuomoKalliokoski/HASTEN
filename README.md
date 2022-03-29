```

    )        (                 )
 ( /(  (     )\ )  *   )    ( /(
 )\()) )\   (()/(` )  /((   )\())
((_)((((_)(  /(_))( )(_))\ ((_)\ 
 _((_)\ _ )\(_)) (_(_()|(_) _((_)
| || (_)_\(_) __||_   _| __| \| |
| __ |/ _ \ \__ \  | | | _|| .` |
|_||_/_/ \_\|___/  |_| |___|_|\_|
```

HASTEN (macHine leArning booSTEd dockiNg) version 1.0, https://github.com/TuomoKalliokoski/HASTEN

Reference: Kalliokoski T. Molecular Informatics 2021, doi:10.1002/minf.202100089

# INTRODUCTION
HASTEN is a tool that makes it easier to run machine learning boosted
virtual screening workflows. It is written in very general Python without
relying in non-standard libraries so it is easy to run in any
Python environment.

HASTEN comes by default configured to use chemprop (free) as the machine 
learning engine and Schrödinger Glide (commercial) as the docking engine. It 
is trivial to replace these, just remember that HASTEN assumes that 
smaller the docking score, the better the pose.

HASTEN is designed to work from mega-databases (millions of molecules) upto 
giga-databases (billion of molecules). Largest database tested so far was
4.1 billion (Enamine REAL). As working with larger databases involves
processing of amounts of data, the process has been split into many smaller
tasks that allow efficient distribution of calculations according to ones
computational resources (either hundreds of cores locally or same amount in 
cloud resources).

# INSTALLATION 
Before you can run HASTEN, make sure that you can use following programs:

Schrödinger Suite:
- Phase (used for conformation generation)
- Glide SP (used for docking)

ChemProp:
- chemprop_train for GPU (for training)
- chemprop_pred for CPU (for predicting)

# HOW TO RUN THE SOFTWARE
You should prepare your database to dock into SMILES file (example, mols.smi). 
It is advisable to generate the SMILES using RDKit (for example, Galaxi
had some broken molecules in it).

You can import molecules into database like this:

```
python3 hasten.py -m mols.db -a import-smiles -s mols.smi
```

NOTE: if database is missing, it will be created. Otherwise molecules will be added to the database.

Plan your screen: critical parameter is the number of molecules docked per 
iteration. For smaller databases (<10M) this should be something like 1% 
(example,-f 0.01). When docking giga-sized databases (>1B) this should be 
something like 0.1%. Pick a name for your screen as well (example, -n vs1). 
In first iteration the compounds are selected randomly from the database 
(-i 1) using 10 CPUs (-c 10) and to be docked with dockingtemplate.in 
(LIGANFDILE will be replaced by HASTEN).

## Iteration 1

```
python3 hasten.py -m mols.db -a pick -i 1 -n vs1 -f 0.01 -c 10 -d dockingtemplate.in
```

This will create SMILES files, .inp files for conformation genration and .in files for docking.  Run the job through Schrödinger Job Control (make sure after each step that all jobs are finished before launching the next one):

```
find vs1_iter1_*.inp -exec nice $SCHRODINGER/pipeline -prog phase_db {} -OVERWRITE -HOST localhost:1 -NJOBS 1 \;
find `pwd` -maxdepth 1 -name "vs1_iter1_*.phdb" -exec echo $SCHRODINGER/phase_database {} export -omae {} \; | sed -e 's:.phdb$: -get 1 -limit 99999999:' | parallel --bar -j 1
find vs1_iter1_*.in -exec nice $SCHRODINGER/glide {} -NICE -OVERWRITE -NJOBS 1 -HOST localhost:1 \;
```

Glide will create docking scores in CSV files and PoseView-maegz files. 
Process the raw docking data from Glide:

```
python3 hasten.py -m mols.db -a import-dock -i 1 -n vs1
```

## Iteration 2 and forward

Prepare input files for ML-training:

```
python3 hasten.py -m mols.db -a train -i 2 -n vs1
```

Train the model using chemprop:

```
chemprop_train --dataset_type regression --target_columns docking_score --data_path train_vs1_iter2.csv --save_dir vs1_iter2 --batch_size 250 --no_cache_mol
```

In order to filter out irrelevant predictions for compounds that will not
be considered in the sorting of predicted score, we need to sample the
database (0.001) to get the mean and standard deviation of predicted
docking scores produced by the previously generated model:

```
python3 hasten.py -m mols.db -a sample-pred -i 2 -n vs1 -c 10
```

To run the predictions on a single multicore machine that has enough RAM in it (note: remember to adjust --checkpoint_dir in later iterations):

```
ls -1 samplepred_vs1_iter2_cpu*.csv | awk '{ print "OMP_NUM_THREADS=1 chemprop_predict --num_workers 0 --no_cache_mol --no_cuda --checkpoint_dir vs1_iter2 --test_path",$1,"--preds_path output_"$1 }' | parallel -j 10 --bar
```

Run the command again, now the program detects that you have output_samplepred_*
files in place:

```
python3 hasten.py -m mols.db -a sample-pred -i 2 -n vs1 -c 10
```

Note "Prediction cutoff:" value. This will be used to filter out predictions.
Cut off is calculated as mean - 0.5 * sd and it is saved to cutoff.txt-file.

Run prediction on all molecules on the database using 10 parallel processes
with 10k chunks:

```    
python3 hasten.py -m mols.db -a pred -i 2 -n vs1 -c 10 -x 10000
```

This will create scripts "pred_vs1_iter2_cpu*.sh" that you can modify to change the calculation GPU/CPU. If you run out of RAM, adjust -x parameter to above.
Loading of the model doesn't take that long.

```
ls -1 pred_vs1_iter2_cpu*.sh | sed -e 's:^:sh :' | parallel -j 10 --bar
```

Import results to the database (just imports everything from PRED* subdirectories):

```
python3 hasten.py -a import-pred -n vs1
```

Pick the best predicted, not yet docked molecules for docking:

```
python3 hasten.py -m mols.db -a pick -i 2 -n vs1 -f 0.01 -c 10 -d dockingtemplate.in
```

After the docking, import the docking scores to the existing hasten_vs1.db:

```
find vs1_iter2_*.inp -exec nice $SCHRODINGER/pipeline -prog phase_db {} -OVERWRITE -HOST localhost:1 -NJOBS 1 \;
find `pwd` -maxdepth 1 -name "vs1_iter2_*.phdb" -exec echo $SCHRODINGER/phase_database {} export -omae {} \; | sed -e 's:.phdb$: -get 1 -limit 99999999:' | xargs -0 bash -c
find vs1_iter2_*.in -exec nice $SCHRODINGER/glide {} -NICE -OVERWRITE -NJOBS 1 -HOST localhost:1 \;
```

Process the raw docking data from Glide:

```
python3 hasten.py -m mols.db -a import-dock -i 2 -n vs1
```

# NOTE: YOU WILL PROBABLY WANT TO DOWNLOAD VERSION 1.1 INSTEAD OF THE VERSION USED IN THE PUBLICATION:
# git clone -b version1.1 https://github.com/TuomoKalliokoski/HASTEN

## Read the version 1.1 README at https://github.com/TuomoKalliokoski/HASTEN/tree/version1.1
HASTEN (macHine leArning booSTEd dockiNg)

Written by Tuomo Kalliokoski <tuomo.kalliokoski at orionpharma.com>

See the article at Molecular Informatics (2021) for additional information (http://doi.org/10.1002/minf.202100089)

***************
* INTRODUCTION
***************
HASTEN is a tool that makes it easier to run machine learning boosted
virtual screening workflows. It is written in very general Python without
relying in non-standard libraries so it is easy to run in any
Python environment

Currently only chemprop is supported as machine learning method, but it is
easy to write Shell-scripts to plug-in your own methods. Glide from
Schrodinger is supported in this version, but the same applies here:
it should be easy to plug-in your own docking program. Do note that the
HASTEN assumes that the smaller the docking score, the better the score.

There is also simulation mode, which allows you to run benchmarks using
existing docking_scores without in reality doing anything in 3D. This mode
can be useful when adjusting the machine learning parameters.

* DESCRIPTION OF THE FILES
1. hasten.py -- main program
2. hasten_import.py -- import data
3. hasten_export.py -- export data
4. hasten_import_simulation.py -- allow simulation data to be used
5. hasten_analyze_simulation.py -- calculate recall on simulated data

6. glide.protocol -- example protocol on how to run Glide
7. simulate.protocol -- example protocol on how to run simulations

8. glide_confgen.py -- glide wrappers
9. glide_confgen.sh
10. glide_docking.py
11. glide_docking.sh

12. simulate_confgen.py -- simulation wrappers
13. simulate_confgen.sh
14. simulate_docking.py
15. simulate_docking.sh

* VERSIONS USED IN THE DEVELOPMENT
- chemprop v1.1.0 (Jan 2020)
- CUDA driver v10.1 and v10.2
- anaconda3, conda 4.9.2
- CentOS 7.6.1810 and 7.8.2003

Tested also Azure cloud VM with Tesla V100 and CUDA 11.1 plus
AWS cloud VM with Tesla V100 and CUDA 10.1.

**********************************
* INSTALLING CHEMPROP ON CENTOS 7
**********************************
NOTE: Please see chemprop webpage for up-to-date instructions. Here is just
what I did to get the program running on Jan 2020.

NOTE2: Even if you have chemprop already installed, check scripts
ml_chemprop_train.sh and ml_chemprop_pred.sh and adjust correct GPU ID for
your calculation card on multiple GPU systems!

0. Check the CUDA version of your system (nvidia-smi).

1. Install anaconda3 to your system.

2. "git glone https://github.com/chemprop/chemprop.git"

3. "cd chemprop"

4. Edit "environment.yml" => change python=3.7.9, add 
   cudatoolkit-<your-cuda-version> and set correct PyTorch version
   (must be newer than 1.5.1).

5. "conda env create -f environment.yml"

6. "conda activate chemprop"

7. "pip install -e ."

8. "pip install git+https://github.com/bp-kelley/descriptastorus"

9. Define your GPU ID number to two HASTEN files: 
    "ml_chemprop_pred.sh" and "ml_chemprop_train.sh" [--gpu]. You may
    check your GPU ID numbers with "nvidia-smi" command. If you have
    several computers, it is good idea to use different copies for each
    computer (you may define this file in protocol file).

************************
* HOW TO RUN SIMULATION
************************
PREPARING

1. You should have the SMILES of the whole database (for example, mols.smi)
and the docking scores for each (for example, dock.txt). The text file should
have the docking_score followed by the docking score (delimiter space).

2. Import simulation data as "python hasten_import_simulation.py -s mols.smi -d dock.txt -o testscreen.db". This will create "testscreen.db".

SCREENING PROTOCOL FILE

This file is just a simple text file (comments start with #). See example
"simulate.protocol". Note that you may have to adjust machine learning
shell scripts in multiple GPU-systems (by default, GPU ID 0 is used for
calculations). In any case, edit the shell scripts so that the paths are
correct!

RUNNING A SIMULATION

python hasten.py -m testscreen.db -p simulate.protocol

EXPORTING SIMULATION RESULTS

After hasten.py finished, you may export the final results by typing:

python hasten_export.py -m testscreen.db -c 10.0 -x scores.txt

ANALYZING SIMULATION RESULTS

To calculate recalls at top 1%:

python hasten_analyze_simulation -m testscreen.db -d dock.txt

******************************
* HOW TO RUN HASTEN WITH GLIDE
******************************

PREPARING

1. You should have the SMILES of the whole database (for example, mols.smi).

    python hasten_import.py -o realscreen.db -s mols.smi -d dock.txt

2. You can also include docking scores in text format (dock.txt) as shown
above, but it is optional.

Do note that you can also split the database into small pieces and then
import them file-by-file (handy if you are importing something like Enamine
REAL).

SCREENING PROTOCOL FILE

See example "glide.protocol". Remember to adjust machine learning shell
scripts on multiple GPU machines and paths. Please also see "example.in". 
Note that it is imporant to have these fields in your .in file (the script
will replace INPUTMAEGZ with the input file when running):

LIGANDFILE  INPUTMAEGZ
POSE_OUTTYPE    ligandlib

RUNNING SCREEN

python hasten.py -m realscreen.db -p glide.protocol

********************
* HAND-OPERATED MODE
********************

When working with large (100M+) databases, even the ML calculations start
to take very long time and several computers are needed.

Hand-operated mode allows you to run the process one step at the end and
in parallel calculations in mind.

TWO ITERATION EXAMPLE

db.db => your HASTEN database
para_simulate.protocol => Your simulation protocol

python hasten.py -m db.db -p para_simulate.protocol --hand-operate dock -i 1
python hasten.py -m db.db -p para_simulate.protocol --hand-operate train -i 2
python hasten.py -m db.db -p para_simulate.protocol --hand-operate split-pred -i 2

Copy PRED1-PRED12 to another computers and iter2 model also

At each computer:
python hasten.py -p para_simulate.protocol --hand-operate pred -i 2

After finished, copy *output* into one directory back where db.db is

python hasten.py -m para_simulate.protocol --hand-operate import-pred
python hasten.py -m para_simulate.protocol --hand-operate dock -i 2

Now you have docked two iterations, continue with training iter3 model

******
* TIPS
******
- you may want to start from some other iteration than 1 sometimes. This 
can be defined by using "-i" parameter for hasten.py

- "smiles_confgen.sh" allows you skip the conformer generation completely when
doing docking (useful when your docking script can take SMILES input directly).

- running long runs distributed across different computers is better done
via the hand-operated mode

***************************************************
* INPUT/OUTPUT DATA FORMATS FOR ADDITIONAL PLUG-INS
***************************************************
Conformer generation:

    - hasten.py starts the confgen-script giving two parameters:
        parameter #1: SMILES-file of the compounds. Each compound name is
                      formatted as "<smilesid>|<hastenid>". SMILES is your
                      own molecule ID and hastenid is integer and used for
                      primary key in all Hasten tables.
        parameter 2: the name of HASTEN-db file.

    The conformer script should directly add conformers to "confs"-table
    for the compounds (see glide_confgen.py as an example). You should store
    all forms of the molecule as one big blob to the confs.

Docking:

    - hasten.py starts the docking-script with three parameters:
        parameter #1: the conformations to be docked
        parameter #2: the name of HASTEN .db-file
        parameter #3: smilesid to hastenid mapping, seperated by |
        parameter #4: iteration (integer)

    See scripts "glide_docking.sh" and "glide_docking.py" on how to import
    both dock_score and pose to correct place.

Machine learning training:

    - hasten.py starts the ML-train-script with three parameters:
        parameter #1: training set CSV
        parameter #2: validation set set CSV
        parameter #3: test set CSV
        parameter #4: iteration as "iter1", "iter2", etc...

    Files all have same format: SMILES, hastenid, dock_score.
    This is simple as you don't have to import anything back.

Machine learning prediction:

    - hasten.py submits chunks of predicted molecules in with following
    command line parameters:
        parameter #1: molecules in SMILES,hastenid-format
        parameter #2: iteration in "iter1","iter2", etc.-format

    - the input it expects back must be comma(,)-delimited file:
            column #1: predicted docking score
            column #2: hastenid

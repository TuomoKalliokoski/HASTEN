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
# HASTEN - wrapper for LigPrep
#
# Path for glide_docking.py
glide_docking_py=/data/tuomo/PROJECTS/HASTEN/glide_docking.py
# Number of parallel tasks
cpu=10
# .in filename
infile=example.in
# failed docking score given to compounds that do not produce poses [10.0]
failed_dock_score=10.0
# docking score cutoff for keeping the poses
store_poses_cutoff=-10.0
#
#
tmpin=`mktemp hasten_glideXXXXXX.in`
logfile=`echo $tmpin | sed -e 's:.in:.log:'`
libfile=`echo $tmpin | sed -e 's:.in:_lib.maegz:'`
subjoblog=`echo $tmpin | sed -e 's:.in:_subjobs.log:'`
subjoblogposes=`echo $tmpin | sed -e 's:.in:_subjob_poses.zip:'`
subjobtar=`echo $tmpin | sed -e 's:.in:_subjobs.tar.gz:'`
# rename out to .mae
maename=`echo $1 | sed -e 's:.out:.mae:'`
mv $1 $maename
sed -e "s:INPUTMAEGZ:$maename:" $infile >$tmpin
$SCHRODINGER/glide -HOST localhost:$cpu -NJOBS $cpu -WAIT -OVERWRITE $tmpin
$SCHRODINGER/run $glide_docking_py $libfile $maename $failed_dock_score $2 $3 $store_poses_cutoff
rm -f $tmpin $logfile $subjoblog $subjoblogposes $subjobtar $libfile $maename

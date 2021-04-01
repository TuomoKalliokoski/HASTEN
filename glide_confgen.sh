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
# Where glide_confgen.py is located
glide_confgen_py=/data/tuomo/PROJECTS/HASTEN/glide_confgen.py
# Number of parallel tasks
cpu=20
# Force field (16 = OPLS3, 14 = OPLS2005)
ff=16
#
#
tmpmaegz=`mktemp hasten_ligprepXXXXXX.maegz`
prefix=`echo $tmpmaegz | cut -d. -f1`
logfile=`echo $tmpmaegz | sed -e 's:.maegz:.log:'`
$SCHRODINGER/ligprep -bff $ff -HOST localhost:$cpu -NJOBS $cpu -WAIT -ismi $1 -omae $tmpmaegz -epik
$SCHRODINGER/run $glide_confgen_py $tmpmaegz $2
rm -f $tmpmaegz $logfile $prefix-dropped.*

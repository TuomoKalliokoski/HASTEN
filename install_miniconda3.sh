#!/bin/bash
#
# just install fresh anaconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
__conda_setup="$('/wrk/tuomo/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/wrk/tuomo/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/wrk/tuomo/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/wrk/tuomo/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
conda install conda-token -n root -y
conda token set MYCODE
conda install -c conda-forge -n base mamba -y

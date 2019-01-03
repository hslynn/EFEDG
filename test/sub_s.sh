#!/usr/bin/bash

module load intel/2018_update1
read -p "Input the name of your task:" name
read -p "Input how many procs you want to run:" np
bsub -J ${name} -R "span[ptile=36]" -n ${np} -e ${name}".err" -o ${name}".out" "mpijob -t mvapich2 ./simplest"

#!/usr/bin/bash

module purge
module load GCC/12.2.0 OpenMPI/4.1.4 BUSCO/5.7.1

ref=$1
thread=$2
lineage=$3
modes=$4
export NUMEXPR_MAX_THREADS=$thread
busco -i $ref -c $thread -l $lineage -m $modes -o $ref.busco_out --offline

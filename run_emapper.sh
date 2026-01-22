#!/usr/bin/bash

module purge
module load GCC/12.3.0  OpenMPI/4.1.5 eggnog-mapper/2.1.12

input=$1
output=$2
thread=$3

emapper.py \
	--cpu $thread \
	-i $input \
	--dbmem \
	--go_evidence all \
	--output $output \
	--pfam_realign realign \
	--data_dir /scratch/data/bio/eggnog-mapper/latest

#rm -r emappertmp_*

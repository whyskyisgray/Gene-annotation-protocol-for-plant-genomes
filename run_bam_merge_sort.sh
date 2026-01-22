#!/usr/bin/bash
module purge
module load GCC/13.2.0 SAMtools/1.21

file_list=$1
threads=$2
output=$3

samtools merge -@ $threads -b $file_list - \
	| samtools sort -@ $threads -o $output - 

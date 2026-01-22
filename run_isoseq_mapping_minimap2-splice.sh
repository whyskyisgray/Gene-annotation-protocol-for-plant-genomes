#!/usr/bin/env bash

set -Eeuo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reference.fasta> <reads.fasta[.gz]> <threads>" >&2
    exit 1
fi

reference="$1"
infile="$2"
threads="$3"

module purge
module load GCC/13.2.0  OpenMPI/4.1.6 pgenomics/1.0.1

### --------------------------------------------------------------------------
### 3.  derive prefix and output paths
### --------------------------------------------------------------------------
# directory in which the input file resides
indir=$(dirname "$infile")

# strip only a terminal .fasta.gz  or  .fasta  extension
filename=$(basename "$infile")
prefix=${filename%.fasta.gz}
prefix=${prefix%.fasta}

outbam="${indir}/${prefix}.sorted.bam"

mm2plus -ax splice:hq -t "$threads" "$reference" "$infile" \
  | samtools sort -@ "$threads" -m 4G -o "$outbam" -
samtools index -@ "$threads" "$outbam"



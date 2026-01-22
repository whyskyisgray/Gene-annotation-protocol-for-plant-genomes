#!/usr/bin/env bash
###############################################################################
# split_into_N_parts.sh  ―  auto-width numeric suffix
#
# Usage :  split_into_N_parts.sh <input_file> <N_parts> [prefix]
# Example: split_into_N_parts.sh reads.fastq 100       # → reads.fastq.part_001 … 100
#          split_into_N_parts.sh reads.fastq 10000 fq_ # → fq_00001 … 10000
###############################################################################

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <input_file> <N_parts> [prefix]" >&2
  exit 1
fi

infile=$1
nparts=$2
prefix=${3:-"$(basename "$infile").part_"}

[[ ! -f $infile ]]                && { echo "Error: '$infile' not found." >&2; exit 2; }
[[ ! $nparts =~ ^[1-9][0-9]*$ ]]  && { echo "Error: N_parts must be > 0." >&2;  exit 3; }

# Width = number of digits in N_parts (100 → 3; 10000 → 5, …)
width=${#nparts}

total=$(wc -l < "$infile")
per=$(( (total + nparts - 1) / nparts ))   # ceiling division

split -l "$per" \
      -d -a"$width" --numeric-suffixes=1 \
      "$infile" "$prefix"

echo "✓ Split into $nparts files with ${width}-digit suffixes (prefix ‘$prefix’)"

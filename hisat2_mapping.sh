#!/usr/bin/env bash

set -euo pipefail

threads="$1"; shift
index="$1";   shift
prefix="$1";  shift
reads_dir="$1"

# ────────────────────────── modules ─────────────────────────────────────────
module purge
module load GCC/13.2.0 OpenMPI/4.1.6 \
            HISAT2/2.2.1 SAMtools/1.21 \
            Trimmomatic/0.39-Java-11

# ───────────────────── scratch in reads_dir ─────────────────────────────────
scratch=$(mktemp -d "${reads_dir}/${prefix}.tmp.XXXXXX")
trap 'rm -rf "$scratch"' EXIT        # always clean up

r1_trim="${scratch}/${prefix}.R1.trimmed.fq.gz"
r2_trim="${scratch}/${prefix}.R2.trimmed.fq.gz"

log_trim="${reads_dir}/${prefix}.trimmomatic.err"
log_hisat="${reads_dir}/${prefix}.hisat2.err"
log_sort="${reads_dir}/${prefix}.samtools_sort.err"
bam_out="${reads_dir}/${prefix}.sorted.bam"

# ─────────────────────── Trimmomatic (zcat stream) ─────────────────────────
trimmomatic PE -threads "$threads" -phred33 \
  <(find "$reads_dir" -regextype posix-extended \
      -regex '.*(_R1|_1|r1).*f(ast)?q\.gz$' -print0 | xargs -0 zcat) \
  <(find "$reads_dir" -regextype posix-extended \
      -regex '.*(_R2|_2|r2).*f(ast)?q\.gz$' -print0 | xargs -0 zcat) \
  "$r1_trim" /dev/null \
  "$r2_trim" /dev/null \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
  2> "$log_trim"

# ─────────────────────── HISAT2 → samtools sort ────────────────────────────
hisat2 -p "$threads" -x "$index" \
       --rna-strandness RF --max-intronlen 10000 \
       --no-mixed --no-discordant \
       -1 "$r1_trim" -2 "$r2_trim" \
       2> "$log_hisat" | \
samtools sort -@ "$threads" -o "$bam_out" - \
       2> "$log_sort"

# ─────────────────────── remove trimmed FASTQ ───────────────────────────────
rm -f "$r1_trim" "$r2_trim"

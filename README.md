# BRAKER based gene annotation pipeline 


---

## 0. Paths and variables

### Code and reference

```bash
codes="/scratch/group/pgenomics/share/temp_codes"
development="/scratch/group/pgenomics/share/development"
ref="/scratch/group/pgenomics/PG1/jayakodi/zoysia/zj_palisades/haphic/haphic_results/final_seq/250411_ZJ_palisades_haplo_resolved_v1.fasta"
```

---

## 1. Short read RNA seq mapping with HISAT2

### Modules

```bash
module purge
module load GCC/13.2.0 OpenMPI/4.1.6 HISAT2/2.2.1
```

### Build index

```bash
mkdir 00_illumina_mapping
index_name="$ref.db"
hisat2-build -p 8 "$ref" "$index_name"

db_name=$(realpath "$index_name")
```

### Submit SRA sample mappings

```bash
cat sample.ids | while read line
do
  echo "sbatch -J hisat -o $line.out -e $line.err -N 1 -c 20 --mem 100G --time 12:00:00 $codes/hisat2_mapping.sh 20 $db_name sra $line"
done | parallel
```

### Submit internal sample mappings

```bash
cat internal_sample.ids | while read line
do
  echo "sbatch -J hisat -o $line.out -e $line.err -N 1 -c 8 --mem 100G --time 1- $codes/hisat2_mapping.sh 8 $db_name sra $line"
done | parallel
```

### Merge and sort BAM

Single step merge and sort example

```bash
sbatch -N 1 -c 30 --mem 150G --time 1- $codes/run_bam_merge_sort.sh bam_files 30 sr.merge.sort.bam
sr_bam=`realpath sr.merge.sort.bam`
```

---

## 2. Iso Seq FASTA preparation and selection of complete ORF transcripts

### Modules

```bash
module purge
module load GCC/13.2.0 SAMtools/1.21
```

### Extract FASTA ids and split

```bash
isoseq_hq_fasta="all.fasta"
prefix="iso"

grep ">" "$isoseq_hq_fasta" | cut -c 2- > fa.ids

bash $pgenomics/share/temp_codes/split_into_N_parts.sh fa.ids 8 "$prefix"

ls -d * | grep "$prefix" | while read line
do
  echo "xargs samtools faidx $isoseq_hq_fasta < $line > $line.fasta"
done | parallel
```

### TransDecoder to identify transcripts with complete ORFs

```bash
module purge
module load GCC/11.3.0 OpenMPI/4.1.4 TransDecoder/5.7.1

ls "$prefix"*.fasta | while read line
do
  echo "TransDecoder.LongOrfs -t $line"
done | parallel
```

Generate complete ORF transcript id lists

```bash
find . -name "*longest_orfs.pep" | while read line
do
  cat "$line" \
    | grep complete \
    | cut -d ' ' -f 1 \
    | cut -c 2- \
    | sed 's/\.p/\t/g'
    | cut -f 1 \
    | sort -u > "$line.complete.transcript.ids"
done
```

Extract complete transcript sequences

```bash
module purge
module load GCC/13.2.0 SAMtools/1.21

find . -name "*.complete.transcript.ids" | while read line
do
  echo "xargs samtools faidx $isoseq_hq_fasta < $line | fold -60 > $line.fasta"
done | parallel -j 8
```

---

## 3. Iso Seq mapping with minimap2 splice

### Submit mappings

```bash
mkdir 01_isoseq_mapping

cat isoseq.ids | while read line
do
  echo "sbatch -e $line.map.err -o $line.map.err -J iso_map -N 1 -c 20 --mem 300G --time 1- $codes/run_isoseq_mapping_minimap2-splice.sh $ref $line 20"
done | parallel
```

### Merge BAM files

```bash
samtools merge -@ 24 -b iso_bam.list - \
  | samtools sort -@ 8 -o iso_merge.sort.bam -

sbatch -N 1 -c 20 --mem 200G --time 1- $codes/run_bam_merge_sort.sh iso_bam.files 20 iso_bam.merge.bam
```

---

## 4. BRAKER run with short read and protein evidence

Recommended to run in tmux

```bash
tmux new -s braker
srun -N 1 -c 30 --mem 200G --time 1- --pty bash

module purge
module load GCC/12.3.0 OpenMPI/4.1.5 BRAKER/3.0.8-long_reads
```

Run BRAKER

```bash
sp_name="zm_sr"

braker.pl \
  --genome="$ref" \
  --species="$sp_name" \
  --bam="$sr_bam" \
  --prot_seq=/scratch/group/pgenomics/PG1/lee/01_turf/03_gene_annotation/zj/annotation_v1/02_allied_prot/all.pep.fa \
  --gff3 \
  --gm_max_intergenic 10000 \
  --AUGUSTUS_CONFIG_PATH=/scratch/group/pgenomics/PG1/lee/pgenomics_biotools/BRAKER/scripts/cfg/ \
  --threads=30
```

---

## 5. Long read model generation with Cupcake, GeneMarkS T, and coordinate conversion

### Isoform collapse

```bash
module purge
module load GCC/11.3.0 OpenMPI/4.1.4 cDNA_Cupcake/29.0.0

long_read_dir="/scratch/group/pgenomics/PG1/lee/01_turf/03_gene_annotation/zj/01_iso-seq_mapping/isoseq_bam_to_fasta/04_hq_fasta"

collapse_isoforms_by_sam.py \
  --input <(zcat $long_read_dir/*fasta.gz) \
  -s <(samtools view iso_merge.sort.bam) \
  --dun-merge-5-shorter \
  -o cupcake
```

### Convert Cupcake GFF to FASTA and run GMST

```bash
module purge
module load GCC/12.3.0 OpenMPI/4.1.5 BRAKER/3.0.8-long_reads

stringtie2fa.py -g "$ref" -f cupcake.collapsed.gff -o cupcake.fa

gmst.pl \
  --strand direct cupcake.fa.mrna \
  --output gmst.out \
  --format GFF

gmst2globalCoords.py \
  -t cupcake.collapsed.gff \
  -p gmst.out \
  -o gmst.global.gtf \
  -g "$ref"
```

---

## 6. Evidence integration with TSEBRA

```bash
module purge
module load GCC/12.3.0 OpenMPI/4.1.5 BRAKER/3.0.8-long_reads

wdir="/scratch/group/pgenomics/PG1/lee/01_turf/03_gene_annotation/zj/annotation_v1/03_annotation"

tsebra.py \
  -g short_read_prot/braker/Augustus/augustus.hints.gtf \
  -e short_read_prot/braker/hintsfile.gff \
  -l gmst.global.longest.gtf \
  -c $working_dir/pgenomics_biotools/TSEBRA/config/long_reads.cfg \
  -o tsebra.gtf -kl
```

---

## 7. GFF3 cleanup, longest isoform, ID normalization, UTR removal, and model merging

### Format normalization using AGAT

```bash
module purge
module load GCC/13.2.0 AGAT/1.4.2

ls gmst.global.gtf tsebra.gtf | rev | cut -d . -f 2- | rev | while read line
do
  agat_sp_keep_longest_isoform.pl -gff $line.gtf -out $line.longest.gff
  agat_convert_sp_gxf2gxf.pl --gff $line.longest.gff -o $line.longest.fix.gff
done

agat_sp_manage_IDs.pl --gff tsebra.longest.fix.gff --prefix tsebra_ --tair -o tsebra.longest.rename.gff3 > /dev/null 2> /dev/null &
agat_sp_manage_IDs.pl --gff gmst.global.longest.fix.gff --prefix long_ --tair -o gmst.global.longest.rename.gff3 > /dev/null 2> /dev/null &
```

Normalize transcript feature names and sources

```bash
cat tsebra.longest.rename.gff3 | sed 's/;gene_id/,/g' | cut -d ',' -f 1 > tsebra.longest.rename.2.gff3
cat gmst.global.longest.rename.gff3 | sed 's/;gene_id/,/g' | cut -d ',' -f 1 > gmst.global.longest.rename.2.gff3

sed -i 's/mRNA/transcript/g' *2.gff3
sed -i 's/AGAT/AUGUSTUS/g' tsebra.longest.rename.2.gff3
sed -i 's/AGAT/GeneMark.hmm/g' gmst.global.longest.rename.2.gff3
```

Resolve overlaps and keep longest isoforms

```bash
agat_sp_fix_overlaping_genes.pl -f tsebra.longest.rename.2.gff3 -o tsebra_fix_overlap.gff3
agat_sp_keep_longest_isoform.pl -f tsebra_fix_overlap.gff3 -o tsebra_fix_overlap.longest.gff3
```

Remove UTR features while touching transcripts

```bash
python $development/gene_annotation/remove_utr_in_gff3.py \
  -i tsebra_fix_overlap.longest.gff3 --touch-transcripts \
  | grep -vi utr > tsebra_no_utr.gff3

python $development/gene_annotation/remove_utr_in_gff3.py \
  -i gmst.global.longest.rename.2.gff3 --touch-transcripts \
  | grep -vi utr > gmst_no_utr.gff3
```

Merge long read and short read evidence

```bash
python $development/gene_annotation/long_short_gff_evidence_merge-sh3.py \
  -l gmst_no_utr.gff3 \
  -s tsebra_no_utr.gff3 \
  -o annot_edit \
  --min_short_overlap 0.7
```

---

## 8. BUSCO evaluation and optional rescue

```bash
module purge
module load GCCcore/12.3.0 gffread/0.12.7

gffread -g "$ref" -y annot_edit_busco_check_prot.pep annot_edit_merged.gff3

sbatch -N 1 -c 48 --mem 100G --time 3:00:00 \
  /scratch/user/saehyun.lee/codes/run_busco.sh \
  annot_edit_busco_check_prot.pep 48 \
  $pgenomics/share/busco_downloads/poales_odb10 prot
```

If BUSCO is low, re evaluate tsebra models

```bash
gffread -g "$ref" -y tsebra_no_utr_busco_check.pep tsebra_no_utr.gff3

sbatch -N 1 -c 48 --mem 100G --time 3:00:00 \
  /scratch/user/saehyun.lee/codes/run_busco.sh \
  tsebra_no_utr_busco_check.pep 48 \
  $pgenomics/share/busco_downloads/poales_odb10 prot
```

Extract BUSCO keep ids

```bash
find annot_edit_busco_check_prot.pep.busco_out/run_poales_odb10/busco_sequences -name "*.faa" \
  | xargs cat \
  | grep ">" \
  | cut -c 2- \
  | rev | cut -d . -f 2- | rev > busco.keep.ids
```

---

## 9. TE overlap based filtering with functional gene rescue

### Identify genes with at least 90 percent overlap with TE annotations

```bash
module purge
module load GCC/13.2.0 BEDTools/2.31.1

awk '$3=="gene"' annot_edit_merged.gff3 > gene.gff
repeat_gff="/scratch/group/pgenomics/PG1/lee/01_turf/04_repeat_masking/ZM/20250710_zm_diamond_hap_resolvd_pseudomolecules_v2.fasta.mod.EDTA.TEanno.gff3"

bedtools intersect -a gene.gff -b "$repeat_gff" -wo \
  | awk '{ gene_len=$5-$4+1; if (($NF/gene_len)>=0.9) print $1"\t"$4"\t"$5 }' \
  | sort -u > repeat_gene_intersect_region.bed
```

### Remove TE overlapped genes using AGAT

```bash
module purge
module load GCC/13.2.0 AGAT/1.4.2

agat_sp_filter_record_by_coordinates.pl \
  --gff annot_edit_merged.gff3 \
  --ranges repeat_gene_intersect_region.bed \
  -o genes.noTE > genes.noTE.log
```

Downstream steps use eggNOG mapper annotations and BUSCO keep ids to rescue functional genes from repeat regions, while removing TE related domains.

---

## 10. lncRNA filtering with FEELnc

Convert to GTF

```bash
module purge
module load GCCcore/12.3.0 gffread/0.12.7

gffread -T -o no_repeat.gtf no_repeat.gff3
```

Run FEELnc

```bash
module purge
ml FEELnc/0.2

mkdir -p lnc_RNA_annot
cp no_repeat* lnc_RNA_annot
cd lnc_RNA_annot

FEELnc_filter.pl -i no_repeat.gtf -a no_repeat.gtf -b transcript_biotyupe=protein_coding --monoex=1 > candidates.gtf

export FEELNCPATH=/sw/eb/sw/FEELnc/0.2/

FEELnc_codpot.pl \
  -i candidates.gtf \
  -a no_repeat.gtf \
  -b transcript_biotype=protein_coding \
  -g "$ref" \
  --mode=shuffle

FEELnc_classifier.pl \
  -i feelnc_codpot_out/candidates.gtf.lncRNA.gtf \
  -a no_repeat.gtf \
  > lncRNA_classes.txt

cut -f 2 lncRNA_classes.txt | awk 'NR>1' | sort -u > lnc.ids
cd ../
```

Remove lncRNA candidates

```bash
module purge
module load GCC/13.2.0 AGAT/1.4.2

agat_sp_filter_feature_from_kill_list.pl \
  -f no_repeat.gff3 -kl lnc_RNA_annot/lnc.ids \
  -o no_lnc_no_repeat.gff3 > agat_lnc_filter.log
```

---

## 11. Remove short proteins and additional TE like models, then build final keep list

This step removes models shorter than 100 aa unless supported by BUSCO, and removes TE related domains from PFAM and annotation fields. A final keep list is generated and used to write the final merged GFF3.

---

## 12. Split haplotypes, BUSCO, renaming, and final outputs

Split by haplotype

```bash
(grep '^#' final_merged.gff3; grep -v '^#' final_merged.gff3 | grep H1 | sort -V -k1,1 -k4,4n) > H1.gff
(grep '^#' final_merged.gff3; grep -v '^#' final_merged.gff3 | grep H2 | sort -V -k1,1 -k4,4n) > H2.gff
```

Extract proteins and CDS and run BUSCO

```bash
module purge
module load GCCcore/12.3.0 gffread/0.12.7

gffread -g "$ref" -y H1.pep H1.gff
gffread -g "$ref" -y H2.pep H2.gff
```

Final naming steps use AGAT manage IDs and a custom renaming script, then final files are produced

- H1.final.gff, H2.final.gff
- H1.final.pep, H2.final.pep
- H1.final.cds, H2.final.cds

---

## 13. One to one ortholog detection between haplotypes by reciprocal best BLAST hits

Inputs

```bash
HAP1_FASTA="H1.final.cds"
HAP2_FASTA="H2.final.cds"
```

Run BLAST and compute reciprocal best hits with identity greater than 95 and same chromosome

```bash
module load GCC/12.2.0 OpenMPI/4.1.4 BLAST+/2.14.0

makeblastdb -in "$HAP1_FASTA" -dbtype nucl -out hap1_db
makeblastdb -in "$HAP2_FASTA" -dbtype nucl -out hap2_db

blastn -query "$HAP1_FASTA" -db hap2_db -out hap1_vs_hap2.txt -num_threads 6 -evalue 1e-5 -outfmt 6
blastn -query "$HAP2_FASTA" -db hap1_db -out hap2_vs_hap1.txt -num_threads 6 -evalue 1e-5 -outfmt 6

awk '!seen[$1]++' hap1_vs_hap2.txt > hap1_best_hits.txt
awk '!seen[$1]++' hap2_vs_hap1.txt > hap2_best_hits.txt

awk 'NR==FNR {a[$1]=$2; next} ($2 in a) && (a[$2] == $1)' hap2_best_hits.txt hap1_best_hits.txt \
  | awk '{split($1,a,"."); split($2,b,"."); if(a[4]==b[4]) print $0}' \
  | awk '$3 > 95' > one_to_one_orthologs.txt
```

---

## 14. Hemizygous gene inference

A gene is treated as hemizygous if it is present in a haplotype CDS list but absent from the one to one ortholog set.

```bash
ls *cds | while read line
do
  grep ">" "$line" | cut -c 2- > "$line.ids"
done

cut -f 1 one_to_one_orthologs.txt > hap1.homologous.ids
cut -f 2 one_to_one_orthologs.txt > hap2.homologous.ids

cat H1.final.cds.ids hap1.homologous.ids | sort | uniq -u | tee H1.hemizygous.ids | wc -l
cat H2.final.cds.ids hap2.homologous.ids | sort | uniq -u | tee H2.hemizygous.ids | wc -l
```

---

## Outputs

- final directory
  - H1.final.gff, H2.final.gff
  - H1.final.pep, H2.final.pep
  - H1.final.cds, H2.final.cds
  - BUSCO outputs, eggNOG mapper outputs

- one_to_one_orthologs.txt
- H1.hemizygous.ids, H2.hemizygous.ids

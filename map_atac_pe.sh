#!/bin/bash

# Added safety. If any command fails or a pipe breaks, the script will stop running.
set -eu -o pipefail -o verbose

R1="$1"
R2="$2"
prefix=$(echo "$R1" | sed 's/.fastq.gz//g')
cores="$3"
index="$4"

# Trim the nextera adapters
sh /home/sbuckberry/working_data_01/bin/bbmap/bbduk2.sh \
in="$R1" \
in2="$R2" \
out="$prefix".tmp_R1.fq \
out2="$prefix".tmp_R2.fq \
rliteral=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG,GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG ktrim=r \
mink=3 \
threads="$cores" \
overwrite=true

## bowtie2 alignment
(/usr/local/packages/bowtie2-2.2.5/bowtie2 -q --threads "$cores" -X2000 \
-x "$index" \
-1 "$prefix".tmp_R1.fq -2 "$prefix".tmp_R2.fq | \
samtools view -bSu - | samtools sort -T sorted - > "$prefix".bam) 2> "$prefix".log

# Create the bam file index
sambamba_v0.5.9 index "$prefix".bam

# Remove the tmp files
rm "$prefix".tmp_R1.fq "$prefix".tmp_R2.fq 



#!/bin/bash
#SBATCH --job-name=fastp_trimming
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mia.nahom@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Trimming/QC of reads using fastp
#################################################################
module load fastp/0.23.2
module load parallel/20180122

# set input/output directory variables
INDIR=/core/cbc/core_projects/2023/Lee/Merged_Samples
REPORTDIR=fastp_reports
mkdir -p $REPORTDIR
TRIMDIR=trimmed_sequences
mkdir -p $TRIMDIR

SAMPLES=/core/cbc/core_projects/2023/Lee/samplelist.txt

# run fastp in parallel, 4 samples at a time
cat $SAMPLES | parallel -j 4 \
fastp \
	--in1 $INDIR/{}_Merged_R1.fastq.gz \
	--in2 $INDIR/{}_Merged_R2.fastq.gz \
	--out1 $TRIMDIR/{}_trim_1.fastq.gz \
	--out2 $TRIMDIR/{}_trim_2.fastq.gz \
	--json $REPORTDIR/{}_fastp.json \
	--html $REPORTDIR/{}_fastp.html

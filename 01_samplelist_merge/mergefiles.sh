#!/bin/bash
#SBATCH --job-name=mergelanes
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=15G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mia.nahom@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date


###### merging files (concatenate)#########

INDIR=/core/cbc/core_projects/2023/Lee/ECL_84tRNA_Bloods_March2023
OUTDIR=/core/cbc/core_projects/2023/Lee/Merged_Samples
mkdir -p $OUTDIR

TMPDIR=/core/cbc/core_projects/2023/Lee/tmp
SAMPLES=$(cat /core/cbc/core_projects/2023/Lee/samplelist.txt)

###### Concatenate files together #########


for i in $SAMPLES
    do echo "Merging R1"

cat $INDIR/"$i"_L002_R1_001.fastq.gz $INDIR/"$i"_L003_R1_001.fastq.gz >> $OUTDIR/"$i"_Merged_R1.fastq.gz

       echo "Merging R2"

cat $INDIR/"$i"_L002_R2_001.fastq.gz $INDIR/"$i"_L003_R2_001.fastq.gz >> $OUTDIR/"$i"_Merged_R2.fastq.gz

done;

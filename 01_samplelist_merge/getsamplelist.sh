#!/bin/bash
#SBATCH --job-name=samplelist
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

find /core/cbc/core_projects/2023/Lee/ECL_84tRNA_Bloods_March2023/ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq >> /core/cbc/core_projects/2023/Lee/samplelist.txt

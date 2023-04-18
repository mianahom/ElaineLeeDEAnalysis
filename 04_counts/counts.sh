#!/bin/bash
#SBATCH --job-name=htseq_count
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 5
#SBATCH --mem=40G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mia.nahom@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Generate Counts 
#################################################################
module load htseq/0.13.5
module load parallel/20180122

INDIR=/core/cbc/core_projects/2023/Lee/03_mapping/alignments
OUTDIR=counts
mkdir -p $OUTDIR

# accession list
SAMPLES=/core/cbc/core_projects/2023/Lee/samplelist.txt

# gtf formatted annotation file
GTF=/isg/shared/databases/alignerIndex/animal/homo_sapiens/current/Homo_sapiens.GRCh38.105.gtf

# run htseq-count on each sample, up to 5 in parallel
cat $SAMPLES | \
parallel -j 5 \
    "htseq-count \
        -s no \
        -r pos \
        -f bam $INDIR/{}.bam \
        $GTF \
        > $OUTDIR/{}.counts"

#!/bin/bash
#SBATCH --job-name=fasta
#SBATCH --output=fasta%j.txt
#SBATCH --error=fasta%j.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=40

bedtools
REGION_BED="/home/twyang2/scenic/Outputs/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="/home/twyang2/ref/mm10.fa"
CHROMSIZES="/home/twyang2/ref/mm10.chrom.sizes"
DATABASE_PREFIX="multiome_1kb_bg_with_mask"
SCRIPT_DIR="/home/twyang2/create_cisTarget_databases"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        mm10.multiome_1kb_bg_with_mask.fa \
        1000 \
        yes
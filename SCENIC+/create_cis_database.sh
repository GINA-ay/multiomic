#!/bin/bash
#SBATCH --job-name=cis_database_
#SBATCH --output=cis_database_%j.txt
#SBATCH --error=cis_database_%j.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=40

OUT_DIR=""${PWD}""
CBDIR="${OUT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
FASTA_FILE="${OUT_DIR}/mm10.multiome_1kb_bg_with_mask.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"
SCRIPT_DIR="/home/twyang2/create_cisTarget_databases"
DATABASE_PREFIX="multiome_1kb_bg_with_mask"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -c ${OUT_DIR}/cbust \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 20

#!/bin/bash

#SBATCH --job-name=cellranger_mapping
#SBATCH --output=cellranger_mapping_%j.txt
#SBATCH --error=cellranger_mapping_%j.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=40

cellranger-arc count --id MCMV_D90 \
--reference /home/twyang2/ref/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
--libraries /home/twyang2/cellranger/counts/MCMV_D90_lib.csv

cellranger-arc count --id MCMV_D7 \
--reference /home/twyang2/ref/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
--libraries /home/twyang2/cellranger/counts/MCMV_D7_lib.csv

cellranger-arc count --id LCMV_C13_D7 \
--reference /home/twyang2/ref/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
--libraries /home/twyang2/cellranger/counts/LCMV_C13_D7_lib.csv

cellranger-arc count --id LCMV_C13_D21 \
--reference /home/twyang2/ref/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
--libraries /home/twyang2/cellranger/counts/LCMV_C13_D21_lib.csv

cellranger-arc count --id LCMV_Arm_D7 \
--reference /home/twyang2/ref/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
--libraries /home/twyang2/cellranger/counts/LCMV_Arm_D7_lib.csv

cellranger-arc count --id LCMV_Arm_D21 \
--reference /home/twyang2/ref/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
--libraries /home/twyang2/cellranger/counts/LCMV_Arm_D21_lib.csv

#!/bin/bash
#SBATCH --job-name=cistopic
#SBATCH --output=cis_%j.txt
#SBATCH --error=cis_%j.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=40
pycistopic qc --fragments /home/twyang2/cellranger/counts/MCMV_D7/outs/atac_fragments.tsv.gz --regions /home/twyang2/scenic/Outputs/consensus_peak_calling/consensus_regions.bed --tss /home/twyang2/scenic/Outputs/qc/tss.bed --output /home/twyang2/scenic/Outputs/qc/MCMV_D7
pycistopic qc --fragments /home/twyang2/cellranger/counts/MCMV_D90/outs/atac_fragments.tsv.gz --regions /home/twyang2/scenic/Outputs/consensus_peak_calling/consensus_regions.bed --tss /home/twyang2/scenic/Outputs/qc/tss.bed --output /home/twyang2/scenic/Outputs/qc/MCMV_D90
pycistopic qc --fragments /home/twyang2/cellranger/counts/LCMV_Arm_D7/outs/atac_fragments.tsv.gz --regions /home/twyang2/scenic/Outputs/consensus_peak_calling/consensus_regions.bed --tss /home/twyang2/scenic/Outputs/qc/tss.bed --output /home/twyang2/scenic/Outputs/qc/LCMV_Arm_D7
pycistopic qc --fragments /home/twyang2/cellranger/counts/LCMV_Arm_D21/outs/atac_fragments.tsv.gz --regions /home/twyang2/scenic/Outputs/consensus_peak_calling/consensus_regions.bed --tss /home/twyang2/scenic/Outputs/qc/tss.bed --output /home/twyang2/scenic/Outputs/qc/LCMV_Arm_D21
pycistopic qc --fragments /home/twyang2/cellranger/counts/LCMV_C13_D7/outs/atac_fragments.tsv.gz --regions /home/twyang2/scenic/Outputs/consensus_peak_calling/consensus_regions.bed --tss /home/twyang2/scenic/Outputs/qc/tss.bed --output /home/twyang2/scenic/Outputs/qc/LCMV_C13_D7
pycistopic qc --fragments /home/twyang2/cellranger/counts/LCMV_C13_D21/outs/atac_fragments.tsv.gz --regions /home/twyang2/scenic/Outputs/consensus_peak_calling/consensus_regions.bed --tss /home/twyang2/scenic/Outputs/qc/tss.bed --output /home/twyang2/scenic/Outputs/qc/LCMV_C13_D21

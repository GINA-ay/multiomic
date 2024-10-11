BATCH --job-name=macs2
#SBATCH --output=peakscalling_%j.txt
#SBATCH --error=peakscalling_%j.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=40

#LCMV_Arm_D21
macs2 callpeak -t /home/twyang2/cellranger/counts/LCMV_Arm_D21/outs/atac_possorted_bam.bam -f BAM -g mm -n LCMV_Arm_D21 --outdir /home/twyang2/macs2/LCMV_Arm_D21
#LCMV_Arm_D7
macs2 callpeak -t /home/twyang2/cellranger/counts/LCMV_Arm_D7/outs/atac_possorted_bam.bam -f BAM -g mm -n LCMV_Arm_D7 --outdir /home/twyang2/macs2/LCMV_Arm_D7
#MCMV_D7
macs2 callpeak -t /home/twyang2/cellranger/counts/MCMV_D7/outs/atac_possorted_bam.bam -f BAM -g mm -n MCMV_D7 --outdir /home/twyang2/macs2/MCMV_D7
#MCMV_D90
macs2 callpeak -t /home/twyang2/cellranger/counts/MCMV_D90/outs/atac_possorted_bam.bam -f BAM -g mm -n MCMV_D90 --outdir /home/twyang2/macs2/MCMV_D90
#LCMV_C13_D21
macs2 callpeak -t /home/twyang2/cellranger/counts/LCMV_C13_D21/outs/atac_possorted_bam.bam -f BAM -g mm -n LCMV_C13_D21 --outdir /home/twyang2/macs2/LCMV_C13_D21
#LCMV_C13_D7
macs2 callpeak -t /home/twyang2/cellranger/counts/LCMV_C13_D7/outs/atac_possorted_bam.bam -f BAM -g mm -n LCMV_C13_D7 --outdir /home/twyang2/macs2/LCMV_C13_D7

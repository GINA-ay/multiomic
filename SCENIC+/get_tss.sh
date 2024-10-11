#!/bin/bash
#SBATCH --job-name=cistopic
#SBATCH --output=cis_%j.txt
#SBATCH --error=cis_%j.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=40

# 检查小鼠基因列表的名称
pycistopic tss gene_annotation_list | grep Mouse

# 创建输出目录
mkdir -p /home/twyang2/scenic/Outputs/qc

# 获取 TSS 信息
pycistopic tss get_tss \
    --output /home/twyang2/scenic/Outputs/qc/tss.bed \
    --name "mmusculus_gene_ensembl" \
    --to-chrom-source ucsc \
    --ucsc mm10

# 查看输出文件的前几行
head /home/twyang2/scenic/Outputs/qc/tss.bed | column -t

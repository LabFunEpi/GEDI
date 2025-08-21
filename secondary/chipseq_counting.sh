#!/bin/bash
#SBATCH -J job1
#SBATCH -o slurm.log
#SBATCH -c 8
#SBATCH -t 01:00:00
#SBATCH --mem 16G
#SBATCH -p med-n16-64g

module load deeptools

cd /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/counting/

multiBamSummary BED-file --BED merged_DE.bed \
  -b FOXA2_ChIP_WT_DE_1.NoDup.bam \
     FOXA2_ChIP_WT_DE_2.NoDup.bam \
     FOXA2_ChIP_TKO_DE_1.NoDup.bam \
     FOXA2_ChIP_TKO_DE_2.NoDup.bam \
  -o results_DE.npz --outRawCounts counts_DE.tsv -p 8

multiBamSummary BED-file --BED merged_PP.bed \
  -b FOXA2_ChIP_WT_PP_1.NoDup.bam \
     FOXA2_ChIP_WT_PP_2.NoDup.bam \
     FOXA2_ChIP_TKO_PP_1.NoDup.bam \
     FOXA2_ChIP_TKO_PP_2.NoDup.bam \
  -o results_PP.npz --outRawCounts counts_PP.tsv -p 8


#!/bin/bash
#SBATCH -J job1
#SBATCH -o slurm.log
#SBATCH -c 8
#SBATCH -t 01:00:00
#SBATCH --mem 16G
#SBATCH -p med-n16-64g

module load deeptools

cd /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/counting/

multiBamSummary BED-file --BED merged.bed \
  -b WT_PP_1.NoDup.bam WT_PP_2.NoDup.bam TKO_PP_1.NoDup.bam TKO_PP_2.NoDup.bam \
  -o results.npz --outRawCounts counts.tsv -p 8

# https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html

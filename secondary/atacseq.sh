
mkdir ~/GEDI/atacseq

cd ~/GEDI/atacseq

gcloud auth login --no-launch-browser

gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/results/atacseq/beds/* .
gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/results/atacseq/bigwigs/* .
gcloud storage cp gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/references/annotations/hg38.rds .

module load bedtools

wc -l *.bed

### Continue analysis in atacseq.R ...

###
###
###

### ... continuing from atacseq.R to create a consensus (MERGE)

bedtools intersect -a WT_PP_1.macs2.bed -b WT_PP_2.macs2.bed -u | wc -l
bedtools intersect -a TKO_PP_1.macs2.bed -b TKO_PP_2.macs2.bed -u | wc -l

bedtools intersect -a WT_PP_1.macs2.bed -b WT_PP_2.macs2.bed -u > WT_PP.macs2.bed
bedtools intersect -a TKO_PP_1.macs2.bed -b TKO_PP_2.macs2.bed -u > TKO_PP.macs2.bed

# Pool WT and TKO peaks together

cat WT_PP.macs2.bed TKO_PP.macs2.bed > pooled.bed
sort -k1,1 -k2,2n pooled.bed > pooled_sorted.bed
bedtools merge -i pooled_sorted.bed > merged.bed

wc -l merged.bed

############################################################################
### Important!!! Run the following in Slurm project SSH (NOT IN RStudio) ###
############################################################################

# Remember to change this folder name to yours
cd /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/

mkdir counting

cd counting

gcloud auth login --no-launch-browser

gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/results/atacseq/bams/* .

cp ~/GEDI/atacseq/merged.bed .

cp ~/GEDI/secondary/atacseq_counting.sh .

sbatch atacseq_counting.sh

squeue --me

# After that completes ... 

cp counts.tsv ~/GEDI/atacseq/

### Continue analysis in atacseq.R ...

###
###
###






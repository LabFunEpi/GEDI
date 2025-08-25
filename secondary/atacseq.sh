
### You'll run this if you've generated results from primary analysis in fslustre with nextflow

# mkdir -p ~/GEDI/demo/atacseq
# cd ~/GEDI/demo/atacseq
# cp -r /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/GEDI/results/beds/* .
# cp -r /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/GEDI/results/bigwigs/* .
# cp -r /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/GEDI/references/annotations/hg38.rds .

### We'll take a shortcut; all the required files from primary analysis are compiled in the demo folder

cd ~/GEDI/demo/atacseq

### Intro to BED files

wc -l *.bed

### Continue analysis in atacseq.R ...

###
###
###

### ... continuing from atacseq.R to create a consensus (MERGE)

module load bedtools

bedtools intersect -a WT_PP_1.macs2.bed -b WT_PP_2.macs2.bed -u | wc -l
bedtools intersect -a TKO_PP_1.macs2.bed -b TKO_PP_2.macs2.bed -u | wc -l

bedtools intersect -a WT_PP_1.macs2.bed -b WT_PP_2.macs2.bed -u > WT_PP.macs2.bed
bedtools intersect -a TKO_PP_1.macs2.bed -b TKO_PP_2.macs2.bed -u > TKO_PP.macs2.bed

# Pool WT and TKO peaks together

cat WT_PP.macs2.bed TKO_PP.macs2.bed > pooled.bed
sort -k1,1 -k2,2n pooled.bed > pooled_sorted.bed
bedtools merge -i pooled_sorted.bed > merged.bed

wc -l merged.bed

### Counting reads per peak in merged.bed

############################################################################
### Important!!! Run the following in Slurm project SSH (NOT IN RStudio) ###
############################################################################

# You can skip this if you want (for the demo); counts.tsv is already in the demo folder

# Remember to change this folder name to yours
cd /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/

mkdir counting

cd counting

gcloud auth login --no-launch-browser

gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/results/atacseq/bams/* .

cp ~/GEDI/demo/atacseq/merged.bed .

cp ~/GEDI/secondary/atacseq_counting.sh .

sbatch atacseq_counting.sh

squeue --me

# After that completes ... 
cp counts.tsv ~/GEDI/demo/atacseq/

### Continue analysis in atacseq.R ...

###
###
###






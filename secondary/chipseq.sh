
### You'll run this if you've generated results from primary analysis in fslustre with nextflow

# mkdir -p ~/GEDI/demo/chipseq
# cd ~/GEDI/demo/chipseq
# cp -r /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/GEDI/results/beds/FOXA2*.macs2.bed .
# cp -r /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/GEDI/results/bigwigs/FOXA2*FC.bw .

### We'll take a shortcut; all the required files from primary analysis are compiled in the demo folder

cd ~/GEDI/demo/chipseq

ls

wc -l FOXA2_ChIP_WT_DE_1.macs2.bed
wc -l FOXA2_ChIP_WT_DE_2.macs2.bed

wc -l FOXA2_ChIP_TKO_DE_1.macs2.bed
wc -l FOXA2_ChIP_TKO_DE_2.macs2.bed

wc -l FOXA2_ChIP_WT_PP_1.macs2.bed
wc -l FOXA2_ChIP_WT_PP_2.macs2.bed

wc -l FOXA2_ChIP_TKO_PP_1.macs2.bed
wc -l FOXA2_ChIP_TKO_PP_2.macs2.bed

module load bedtools

bedtools intersect -a FOXA2_ChIP_WT_DE_1.macs2.bed -b FOXA2_ChIP_WT_DE_2.macs2.bed -u > FOXA2_WT_DE.bed
bedtools intersect -a FOXA2_ChIP_TKO_DE_1.macs2.bed -b FOXA2_ChIP_TKO_DE_2.macs2.bed -u > FOXA2_TKO_DE.bed
bedtools intersect -a FOXA2_ChIP_WT_PP_1.macs2.bed -b FOXA2_ChIP_WT_PP_2.macs2.bed -u > FOXA2_WT_PP.bed
bedtools intersect -a FOXA2_ChIP_TKO_PP_1.macs2.bed -b FOXA2_ChIP_TKO_PP_2.macs2.bed -u > FOXA2_TKO_PP.bed

wc -l FOXA2_WT_DE.bed
wc -l FOXA2_TKO_DE.bed
wc -l FOXA2_WT_PP.bed
wc -l FOXA2_TKO_PP.bed

cat FOXA2_WT_DE.bed FOXA2_TKO_DE.bed | sort -k1,1 -k2,2n > concatenated_sorted.bed
bedtools merge -i concatenated_sorted.bed > merged_DE.bed

cat FOXA2_WT_PP.bed FOXA2_TKO_PP.bed | sort -k1,1 -k2,2n > concatenated_sorted.bed
bedtools merge -i concatenated_sorted.bed > merged_PP.bed

wc -l merged_DE.bed
wc -l merged_PP.bed

### Counting reads per peak

############################################################################
### Important!!! Run the following in Slurm project SSH (NOT IN RStudio) ###
############################################################################

# You can skip this if you want (for the demo); counts files are already in the demo folder

# Remember to change this folder name to yours
cd /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/

mkdir counting

cd counting

gcloud auth login --no-launch-browser

gcloud config set project ml-fpt-rsa-maia-s-p-1367

gcloud storage cp gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/results/chipseq/bams/FOXA2_ChIP* .

cp ~/GEDI/demo/chipseq/merged_DE.bed .
cp ~/GEDI/demo/chipseq/merged_PP.bed .

cp ~/GEDI/secondary/chipseq_counting.sh .

sbatch chipseq_counting.sh

squeue --me

cp counts_DE.tsv ~/GEDI/demo/chipseq/

cp counts_PP.tsv ~/GEDI/demo/chipseq/

### Continue analysis in chipseq.R ...

###
###
###

### ... continuing from chipseq.R - Make heatmaps using DeepTools

module load deeptools

computeMatrix reference-point -S FOXA2_ChIP_WT_DE_FC.bw \
                                 FOXA2_ChIP_WT_GT_FC.bw \
                                 FOXA2_ChIP_WT_PP_FC.bw \
                                 FOXA2_ChIP_TKO_PP_FC.bw \
                              -R FOXA2_down_in_PP.bed \
                              --referencePoint center \
                              -b 5000 -a 5000 -p 16 \
                              --outFileNameMatrix FOXA2 \
                              --outFileName FOXA2.tab.gz

plotHeatmap -m FOXA2.tab.gz \
            --colorList 'white,red' \
            --sortRegions keep \
            --kmeans 3 \
            -out FOXA2.pdf

###

computeMatrix reference-point -S FOXA2_ChIP_WT_DE_FC.bw \
                                 FOXA2_ChIP_WT_GT_FC.bw \
                                 FOXA2_ChIP_WT_PP_FC.bw \
                              -R DE_to_PP.bed \
                                 GTh_to_PP.bed \
                                 GT_to_PPh.bed \
                                 PP_specific.bed \
                              --referencePoint center \
                              -b 3000 -a 3000 -p 16 \
                              --outFileNameMatrix FOXA2_hDHMR \
                              --outFileName FOXA2_hDHMR.tab.gz

plotHeatmap -m FOXA2_hDHMR.tab.gz \
            --colorList 'white,blue' \
            --sortRegions keep \
            -out FOXA2_hDHMR.pdf

###

computeMatrix reference-point -S FOXA2_ChIP_WT_PP_FC.bw \
                              -R ~/GEDI/demo/atacseq/merged.bed \
                              --referencePoint center \
                              -b 5000 -a 5000 -p 16 \
                              --outFileNameMatrix FOXK2_ATAC_all \
                              --outFileName FOXK2_ATAC_all.tab.gz &

computeMatrix reference-point -S FOXA2_ChIP_WT_PP_FC.bw \
                              -R ~/GEDI/demo/atacseq/ATAC_down_in_PP.bed \
                              --referencePoint center \
                              -b 5000 -a 5000 -p 16 \
                              --outFileNameMatrix FOXK2_ATAC_down \
                              --outFileName FOXK2_ATAC_down.tab.gz &

# Run after completion of the previous commands
plotHeatmap -m FOXK2_ATAC_all.tab.gz \
            --colorList 'white,red' \
            --sortRegions descend \
            -out FOXK2_ATAC_all.pdf &

plotHeatmap -m FOXK2_ATAC_down.tab.gz \
            --colorList 'white,red' \
            --sortRegions descend \
            -out FOXK2_ATAC_down.pdf &

### Motif analysis
# http://homer.ucsd.edu/homer/ngs/peakMotifs.html

awk -v OFS="\t" '{print $0,"+"}' ~/GEDI/demo/atacseq/ATAC_down_in_PP.bed > ATAC_down_in_PP_4col.bed

module load homer

findMotifsGenome.pl ATAC_down_in_PP_4col.bed \
  /usr/local/biotools/homer/4.11/data/genomes/hg38 \
  ./motifAnalysis \
  -size 200 \
  -preparsedDir ./preparsed



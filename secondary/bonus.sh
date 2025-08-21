# Combine replicates (bigWig files)

module load wiggletools
module load ucscbigtools

wiggletools write WT_PP.wig mean WT_PP_1.bw WT_PP_2.bw &
wiggletools write TKO_PP.wig mean TKO_PP_1.bw TKO_PP_2.bw &

wigToBigWig WT_PP.wig hg38 WT_PP.bw &
wigToBigWig TKO_PP.wig hg38 TKO_PP.bw &

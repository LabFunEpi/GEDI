### Setup for REGS6700: Genomic and Epigenomic Data Integration ###

#################### Prerequisites (Get access) ###################

# 1. MAIA Lab cloud storage buckets on GCP
# 2. HPC GCP Slurm (with folder in the Lustre filesystem)
# 3. RStudio on HPC GCP

###################################################################

### Remember to open links in Incognito window, 
### so that multiple Google accounts don't conflict!!!

# GCP storage in Browser - https://console.cloud.google.com
# SSH via Browser - https://console.cloud.google.com
# RStudio via Browser - https://rstudio-hpc-01.mayo.edu



### Copy required files (code and references) from MAIA lab bucket to fslustre space

# Remember to change this to your folder in fslustre
cd /fslustre/labs/ext_mohammedismail_wazim_mayo_ed 

gcloud auth login --no-launch-browser
gcloud config set project ml-fpt-rsa-maia-s-p-1367

gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/code/GEDI .
# (or) 
# git clone https://github.com/LabFunEpi/GEDI.git 

gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/references ./GEDI



### Copy required files (code and demo data) from MAIA lab bucket to your home

cd ~

gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/code/GEDI .

gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/demo ./GEDI



######################## Nextflow ############################
### Make sure that you've already done the previous steps! ###
##############################################################

cd /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/GEDI

# Prepare the input data

gcloud auth login --no-launch-browser
gcloud config set project ml-fpt-rsa-maia-s-p-1367

gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/data_subset/rnaseq .

mv rnaseq data

# Run rnaseq Nextflow pipeline
module load nextflow
nextflow run rnaseq.nf -profile local

# Check-out the results
ls results
ls results/counts
less results/counts/Sample1_htseq.tab

# Troubleshooting Nextflow runs
nextflow log

# Running Nextflow pipeline with Slurm (and in background)
nextflow -bg -q run rnaseq.nf

# Checking the Slurm jobs running in the background
squeue --me

# Cleaning up after completion of a pipeline
# First, copy the results to storage bucket for use in secondary analysis
gcloud storage cp -r results/ gs://... 

# Remove everything related to this run
rm -rf .nextflow* work/ results/


### Continue secondary analysis in RStudio

# https://rstudio-hpc-01.mayo.edu/
# https://rstudio-hpc-02.mayo.edu/




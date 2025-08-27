
### You'll run this if you've generated results from primary analysis in fslustre with nextflow

# mkdir -p ~/GEDI/demo/dnamethyl
# cd ~/GEDI/demo/dnamethyl
# cp -r /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/GEDI/results/methylKit/* .
# cp -r /fslustre/labs/ext_mohammedismail_wazim_mayo_ed/GEDI/results/bigwigs/* .

### We'll take a shortcut; all the required files from primary analysis are compiled in the demo folder

gcloud auth login --no-launch-browser
gcloud config set project ml-fpt-rsa-maia-s-p-1367

cd ~
gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/demo ./GEDI

cd ~/GEDI/demo/dnamethyl

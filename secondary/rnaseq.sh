
mkdir ~/GEDI/rnaseq

cd ~/GEDI/rnaseq

gcloud auth login --no-launch-browser

gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/results/rnaseq/counts/*_htseq.tab .


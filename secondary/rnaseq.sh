
mkdir ~/GEDI/demo/rnaseq

cd ~/GEDI/demo/rnaseq

gcloud auth login --no-launch-browser

gcloud config set project ml-fpt-rsa-maia-s-p-1367

gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/results/rnaseq/counts/*_htseq.tab .


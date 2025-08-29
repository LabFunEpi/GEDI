
gcloud auth login --no-launch-browser
gcloud config set project ml-fpt-rsa-maia-s-p-1367

cd ~/GEDI/demo/scrnaseq

gcloud storage cp -r gs://ml-phi-proj-rsa-us-central1-p-15fe/REGS6700/demo/scrnaseq/filtered_gene_bc_matrices .


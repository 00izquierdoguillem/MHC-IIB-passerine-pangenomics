#!/bin/bash
#SBATCH -p shared,edwards
#SBATCH -c 2
#SBATCH -t 1-00:00
#SBATCH --mem=64000
#SBATCH -o logs/0_get_canary_isoforms.out
#SBATCH -e logs/0_get_canary_isoforms.err
#SBATCH --mail-type=END

cd /n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/isoforms_Ser_can

# --- Project variables ---

REF_ACC="GCA_022539315.2"
REF_SHORT="Sercan"
REF_DIR="/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/isoforms_Ser_can"

# --- Download NCBI Chicken Reference ---
module load python/3.10.9-fasrc01
conda activate myenv

cd "$REF_DIR"
if [ ! -f "${REF_DIR}/ncbi_dataset.zip" ]; then
    datasets download genome accession $REF_ACC --include gff3,rna,cds,protein,genome,seq-report,gtf
    unzip ncbi_dataset.zip
fi

REF_NCBI_DATA="${REF_DIR}/ncbi_dataset/data/${REF_ACC}"
REF_GTF=$(ls $REF_NCBI_DATA/*.gtf | head -n1)

awk '$3=="transcript" && /transcript_id/ && /gene_id/' $REF_GTF \
| awk '{match($0, /gene_id "([^"]+)"/, g); match($0, /transcript_id "([^"]+)"/, t); if(length(g) && length(t)) print g[1] "\t" t[1]}' | sort | uniq > "$REF_DIR/isoforms.tsv"
sed -i '1igeneID\ttranscriptID' "$REF_DIR/isoforms.tsv"
ISOFORM_FILE="$REF_DIR/isoforms_Sercan.tsv"

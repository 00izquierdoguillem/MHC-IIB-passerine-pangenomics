#!/bin/bash
#SBATCH -p edwards,shared # partition to submit to 
#SBATCH -c 1      # cores requested 
#SBATCH -t 1-00:00  # time requested in days-hours:minutes
#SBATCH -o logs/2_filter_miniprot_based_blast.out # name of output file
#SBATCH -e logs/2_filter_miniprot_based_blast.err  # name of error file
#SBATCH --mem=1000 #Mb
#SBATCH --mail-type=END # email me whe nmy job is done running
#SBATCH --array=1-20%10

module load python/3.10.9-fasrc01
conda activate myenv

WORKDIR=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/BLAST_MINIPROT_overlap

cd $WORKDIR

# 2. After this job, run script 6 to get the final miniprot annotation of the region

# This script should subset every ind's PAF files to get those annotated in the genes that we have previoiusly delimited

# Furthermore, when coordinates overlap between gene annotations, keep just the 1st one (no need for priority as all genes should have equal importance).

#IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scrub_jay.list)
#IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" haeMex.list)
IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" hemMar.list)

PYSCRIPT=/n/home11/gizquierdoaranega/scripts/2_MHC_project/subset_pafs_and_filter.py
GENES=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/fastOMA/fastOMA.genes_HOGs.Jan26.csv
PAF_SJ=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/pangene/paf
PAF_HF=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/BLAST_MINIPROT_overlap
PAF_HEMI=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/BLAST_MINIPROT_overlap
OUTDIR=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/MHC_pafs

conda run -n myenv python $PYSCRIPT $GENES $PAF_HEMI/$IND.unfiltered.paf $OUTDIR/$IND.MHC.paf

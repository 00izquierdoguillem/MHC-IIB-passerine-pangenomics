#!/bin/bash
#SBATCH -p edwards,shared # partition to submit to 
#SBATCH -c 1      # cores requested 
#SBATCH -t 0-01:00  # time requested in days-hours:minutes
#SBATCH -o logs/4_get_fasta_exon2_%a.HF.out # name of output file
#SBATCH -e logs/4_get_fasta_exon2_%a.HF.err  # name of error file
#SBATCH --mem=1000 #Mb
#SBATCH --mail-type=END # email me whe nmy job is done running 
#SBATCH --array=1-50%10

module load python/3.10.9-fasrc01
conda activate myenv

WORKDIR=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/fasta_sequences/Jan26
BLASTDIR=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/BLAST_MINIPROT_overlap

# Script to create a fasta alignment of all exon2 sequences
cd $WORKDIR

IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" haeMex.list)
blastind=$(sed -n "${SLURM_ARRAY_TASK_ID}p" haeMex.blast.list) 

# 1. Quick thing before we start, make sure that all blast hits are below 1e-10

BLASTDIR=$WORKDIR

cd $BLASTDIR
BLASTFILE=$BLASTDIR/hofi_blastn.Jan_26.length_overlap_filtered.out

#awk '$11 <= 1e-10' $BLASTFILE > prova
#mv prova $BLASTFILE

# 2. First create individual files of the blasthits

cd $WORKDIR

less $BLASTFILE | grep "exon3" | cut -f 4,8,9 | grep $blastind | sed 's/\.1//g' | awk '{ $1 = $1 ".1"; print }' OFS='\t' > $IND.ex3.bed

less $IND.ex3.bed | sort | uniq > prova.$IND
mv prova.$IND $IND.ex3.bed

# 3. Then run seqkit to extract the sequences

ASSDIR=/n/netscratch/edwards_lab/Lab/kelsielopez/MHC_test_scott/hofi/renamed_fastas

cp $ASSDIR/$IND\_genomic_renamed.fa $WORKDIR/assemblies/$IND.p_ctg.fa

seqkit subseq --bed $IND.ex3.bed $WORKDIR/assemblies/$IND.p_ctg.fa > $IND.ex3.fa

cat $IND.ex3.fa >> haeMex.ex3.fa


#!/bin/bash
#SBATCH -p edwards,shared # partition to submit to 
#SBATCH -c 1      # cores requested 
#SBATCH -t 0-01:00  # time requested in days-hours:minutes
#SBATCH -o logs/4_get_fasta_exon2_%a.out # name of output file
#SBATCH -e logs/4_get_fasta_exon2_%a.err  # name of error file
#SBATCH --mem=1000 #Mb
#SBATCH --mail-type=END # email me whe nmy job is done running 
#SBATCH --array=1-25%10

module load python/3.10.9-fasrc01
conda activate myenv

WORKDIR=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/fasta_sequences/Jan26
BLASTDIR=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/BLAST_MINIPROT_overlap

# Script to create a fasta alignment of all exon2 sequences
cd $WORKDIR

IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" hemMar.list)
blastind=$(sed -n "${SLURM_ARRAY_TASK_ID}p" hemMar.blast.list) 

# 1. Quick thing before we start, make sure that all blast hits are below 1e-10

cd $BLASTDIR

#awk '$11 <= 1e-10' all_haplotypes_mhc_classII_5_probes_blastn.length_overlap_filtered.out > prova
#mv prova all_haplotypes_mhc_classII_5_probes_blastn.length_overlap_filtered.out

# 2. First create individual files of the blasthits

BLASTDIR=$WORKDIR

cd $WORKDIR

less $BLASTDIR/hemitriccus_blastn.Jan26.length_overlap_filtered.out | grep "exon2" | cut -f 4,8,9 | grep $blastind > $IND.ex2.bed

# 3. Then run seqkit to extract the sequences

ASSDIR=/n/netscratch/edwards_lab/Everyone/kelsie/mhc/MHC_test_scott

cp $ASSDIR/$IND.p_ctg.fa $WORKDIR/assemblies/

seqkit subseq --bed $IND.ex2.bed $WORKDIR/assemblies/$IND.p_ctg.fa > $IND.ex2.fa

cat $IND.ex2.fa >> hemMar.ex2.fa


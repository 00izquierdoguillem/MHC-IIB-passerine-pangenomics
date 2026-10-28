#!/bin/bash
#SBATCH -p test # partition to submit to 
#SBATCH -c 1      # cores requested 
#SBATCH -t 0-01:00  # time requested in days-hours:minutes
#SBATCH -o logs/3_blast_miniprot_overlap.hemi.out # name of output file
#SBATCH -e logs/3_blast_miniprot_overlap.hemi.err  # name of error file
#SBATCH --mem=1000 #Mb
#SBATCH --mail-type=END # email me whe nmy job is done running

WORKDIR=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/BLAST_MINIPROT_overlap/Jan26
BLASTFILE=$WORKDIR/hemitriccus_blastn.Jan26.length_overlap_filtered.out

cd $WORKDIR

# 0. Get only results with eval above 1e-10

awk '$11 < 1e-10' $BLASTFILE > prova
mv prova $BLASTFILE

# 1. Create column individual 

awk -F'\t' '{
    split($4, a, "#");   # a[1] = prefix, a[2] = hap1 or hap2
    hap = (a[2] == "1") ? ".hap1" : (a[2] == "2") ? ".hap2" : "";
    newcol = a[1] hap;
    print $0 "\t" newcol;
}' $BLASTFILE > $BLASTFILE.renamed

# 2. Divide files per individual


for i in $(seq 1 $(wc -l < $WORKDIR/hemMar.list));
do
ind=$(sed -n "${i}p" $WORKDIR/hemMar.list)
grep "$ind" $BLASTFILE.renamed > $WORKDIR/$ind.blast.out
done

# 3. Since the PAF files for hemMar are filtered, we have to resort to GFF files. First, we'll create files that only keep the "PAF" lines from GFFs (instead of the CDS and mRNA ones)

PAFDIR=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/BLAST_MINIPROT_overlap
#PAFDIR=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/def_pafs

#/n/netscratch/edwards_lab/Lab/kelsielopez/blast_validate_cnv/Oct_2025_fix
#for i in $(seq 1 $(wc -l < $WORKDIR/hemMar.list)); 
#do
#ind=$(sed -n "${i}p" $WORKDIR/hemMar.list)
#grep "#PAF" $PAFDIR/$ind\_aln.gff | cut -f 2- > $WORKDIR/$ind.unfiltered.paf
#done


# 4. Find overlap with miniprot GFF files

for i in $(seq 1 $(wc -l < $WORKDIR/hemMar.list)); 
do
ind=$(sed -n "${i}p" $WORKDIR/hemMar.list)

awk -F'\t' '
BEGIN { OFS = FS }

# ---------- Pass 1: load reference intervals (file2) ----------
FNR==NR {
    split($1, g, ":")
    gene = g[1]

    contig = $6
    start  = $8
    end    = $9

    n[contig]++
    s[contig, n[contig]] = start
    e[contig, n[contig]] = end
    gname[contig, n[contig]] = gene
    next
}

# ---------- Pass 2: process file1 ----------
{
    contig1 = $4
    start1  = $8
    end1    = $9
    hit = "NA"

    # normalize coordinates     
    if (start1 > end1) { tmp=start1; start1=end1; end1=tmp }
    if (start2 > end2) { tmp=start2; start2=end2; end2=tmp }

    for (i = 1; i <= n[contig1]; i++) {
        if (start1 <= e[contig1, i] && s[contig1, i] <= end1) {
            hit = gname[contig1, i]
            print $0, hit
            continue
        }
    }

print $0, hit

}' $PAFDIR/$ind.unfiltered.paf $WORKDIR/$ind.blast.out | sort | uniq | sort -k4,4 -nk 8,8 > $WORKDIR/$ind.blast_miniprot_overlap.Jan26.tsv

done

# By keeping the unique rows, we make sure that we don't have any repeated annotation per blast hits (that is, we can have multiple genes per blast hit, but not multiple hits with the same gene annotation)
# The last sort command is to get results ordered by contig and position

# 4. Add column that says whether each miniprot-annotated gene had been previously identified by fastOMA or not

PYSCRIPT=/n/home11/gizquierdoaranega/scripts/2_MHC_project/filter_blast_overlap.py
GENES=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/fastOMA/fastOMA.genes_HOGs.Jan26.csv

for i in $(seq 1 $(wc -l < $WORKDIR/hemMar.list)); 
do
ind=$(sed -n "${i}p" $WORKDIR/hemMar.list)
conda run -n myenv python $PYSCRIPT $WORKDIR/$ind.blast_miniprot_overlap.Jan26.tsv $GENES $WORKDIR/$ind.blast_miniprot_overlap.temp.tsv $WORKDIR/$ind.blast_miniprot_overlap.filtered.Jan26.tsv
mv $WORKDIR/$ind.blast_miniprot_overlap.temp.tsv $WORKDIR/$ind.blast_miniprot_overlap.Jan26.tsv
done


# 4. Now, if you want to get a summary of all miniprot-annotated genes matching with the hits:

cat *.blast_miniprot_overl* | cut -f 13 | sort | uniq

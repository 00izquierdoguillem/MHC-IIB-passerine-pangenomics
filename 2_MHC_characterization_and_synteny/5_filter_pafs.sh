WORKDIR=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/MHC_pafs
cd $WORKDIR

# 1. Substitute gene_names per ortholog names

for infile in *.paf; do
awk -F'\t' -v OFS='\t' 'NR==FNR {
  map[$1] = $2
  next
}{
  if ($21 in map)
    $21 = map[$21]
  print
}' genes.fastOMA.dict $infile | grep -v "*" > prova
mv prova $infile
done

# 2. Create column that say whether each paf annotation is found in the BLAST hits

BLASTFILE=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/BLAST_MINIPROT_overlap/Jan26/allspecies.BLAST_MINIPROT_OVERLAP.Jan26.tsv

for paffile in *.paf; do

awk -F'\t' '
BEGIN { OFS = FS }

# ---------- Pass 1: load reference intervals (file2) ----------
FNR==NR {
    contig = $4
    start  = $8
    end    = $9

    n[contig]++
    s[contig, n[contig]] = start
    e[contig, n[contig]] = end
    next
}

# ---------- Pass 2: process file1 ----------
{
    contig1 = $6
    start1  = $8
    end1    = $9
    hit = "No_blast_overlap"

    # normalize coordinates     
    if (start1 > end1) { tmp=start1; start1=end1; end1=tmp }
    if (start2 > end2) { tmp=start2; start2=end2; end2=tmp }

    for (i = 1; i <= n[contig1]; i++) {
        if (start1 <= e[contig1, i] && s[contig1, i] <= end1) {
            hit = "Blast_overlap"
            print $0, hit
            continue
        }
    }
   print $0, hit
}
' $BLASTFILE $paffile | sort | uniq | sort -k6,6 -nk 8,8 > prova
mv prova $paffile
done


# 3. Join all paf files and add an ind column

cat *.paf > allinds.miniprot_blast_overlap.paf


awk -F'\t' -v OFS='\t' '
{
    # Case 1: column 6 starts with AC, AI, AW, AA, or CY
    if ($6 ~ /^(AC|AI|AW|AA|CY)/) {
        split($6, a, "#")     # a[1]=prefix, a[2]=hap1/hap2
        hap = (a[2] == "1") ? ".1" :
              (a[2] == "2") ? ".2" : ""
        newcol = a[1] hap
    }

    # Case 2: column 6 starts with HMRG
    else if ($6 ~ /^HMRG/) {
        split($6, a, "#")     # a[1]=prefix, a[2]=1/2
        hap = (a[2] == "1") ? ".hap1" :
              (a[2] == "2") ? ".hap2" : ""
        newcol = a[1] hap
    }

    # Case 3: all other cases
    else {
        newcol = $6
        sub(/_[^_]*\.[^.]*$/, "", newcol)
    }

    print $0, newcol
}' allinds.miniprot_blast_overlap.paf > prova

mv prova allinds.miniprot_blast_overlap.paf

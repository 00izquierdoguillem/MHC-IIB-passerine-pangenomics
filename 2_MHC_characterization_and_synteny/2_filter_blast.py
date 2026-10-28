import pandas as pd
import re
import sys

# 1. First read the blast output file

blast = pd.read_csv("scrubjay_blastn.Jan_26.out", sep="\t", header=None)

# We need to a) in overlaps, keep the best hit by evalue and b) eliminate those blast hits in MHCII2 > 230bp and in MHCII3>250bp

# For purpose a) - keep in mind that in case of overlap, we are going to prioritize exon2 over exon3:

#contig - col3
#start - col7
#end - col8
#eval - col10

# Normalize start and end
blast[11]=blast[[7, 8]].min(axis=1)
blast[12]=blast[[7, 8]].max(axis=1)

# Sort by contig and start
blast = blast.sort_values([3, 11], ascending=[True, True])

# Create the "kept intervals"
blast_filtered = []

# And now run by row
for idx, row in blast.iterrows():
        overlap_found = False           # Default
        for i,k in enumerate(blast_filtered):
                if k[3] == row[3] and not (
                        row[12] <= k[11] or 
                        row[11] >= k[12]
                ):
                        # Overlap with a previously kept interval
                        overlap_found = True
                        if row[0]!=k[0] and row[0]!="Mhc_U23958.1_exon3.fa":
                                blast_filtered[i]=row
                                continue
        if not overlap_found:
                blast_filtered.append(row)

blast_filtered = pd.DataFrame(blast_filtered)

blast=blast_filtered

# For purpose b)

#First initialize the output dataframe as well:
blast_filtered_2=[]

# Using continue, no need for an "else" statement in the end

for idx, row in blast.iterrows():
        if row[0] == "Mhc_U23958.1_exon2.fa":
                if int(row[5])>240:
                        blast_filtered_2.append(row)
                else:
                        continue
        elif row[0] == "Mhc_U23958.1_exon3.fa":
                if int(row[5])>250:
                        blast_filtered_2.append(row)
                else:
                        continue
        else:
                blast_filtered_2.append(row)

blast_filtered = pd.DataFrame(blast_filtered_2)

blast=blast_filtered

# c) Finally remove those hits with the first half of exon 3 (col1 <150)

blast_filtered_3=[]

for idx, row in blast.iterrows():
        if row[0] == "Mhc_U23958.1_exon3.fa":
                if int(row[2])>400:
                        blast_filtered_3.append(row)
                else:
                        blast_filtered_3.append(row)
        else:
                blast_filtered_3.append(row)

blast_filtered = pd.DataFrame(blast_filtered_3)

blast_filtered.to_csv("scrubjay_blastn.Jan26.length_overlap_filtered.out", sep="\t", index=False, header=None)

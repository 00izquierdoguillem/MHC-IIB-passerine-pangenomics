import pandas as pd
import re
import sys

# 1. First read the blast output file

blast = pd.read_csv("hemitriccus_blastn.Jan26.out", sep="\t", header=None)

# We need to a) in overlaps, keep the best hit by evalue and b) eliminate those blast hits in MHCII2 > 230bp and in MHCII3>250bp

# For purpose a)

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
        for k in blast_filtered:
                if k[3] == row[3] and not (
                        row[12] <= k[11] or 
                        row[11] >= k[12]
                ):
                        # Overlap with a previously kept interval
                        overlap_found = True
                        # Keep the one already there (higher score due to sorting)
                        break
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
                blast_filtered_2.append(row)

blast_filtered = pd.DataFrame(blast_filtered_2)

blast_filtered.to_csv("hemitriccus_blastn.Jan26.length_overlap_filtered.out", sep="\t", index=False, header=None)

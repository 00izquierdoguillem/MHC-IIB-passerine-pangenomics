import pandas as pd
import re
import sys

# Script to subset the paf file to get only the MHC genes

if len(sys.argv) != 4:
    print("Usage: python subset_pafs.py <gene_list> <input.paf> <output.paf>")
    sys.exit(1)

ortholog_list = sys.argv[1]
input_paf = sys.argv[2]
output_paf = sys.argv[3]

# 0. Load ortholog list file

gene_table = pd.read_csv(ortholog_list, sep="\t")

print("Load ortholog completed")

# 1. Load paf file

paf_hap1 = pd.read_csv(input_paf, sep="\t", header=None, on_bad_lines="warn")

print("Load paf completed")

# Create paf row that holds gene name:

paf_hap1[20] = paf_hap1[0].str.extract(r'([^:]+):')

print("Name column created")

# 2. Get those genes that match and output - this is not doing partial matches. Revise whether this might be a problem for chiLan or wether this could have been a problem in previous steps

paf_hap1_MHC = paf_hap1[paf_hap1.iloc[:,20].isin(gene_table.iloc[:,0])]

print ("Match gene name done")

# 3. Now, for every line, recover only those that don't overlap with other gene annotations

# Normalize start and end
paf_hap1_MHC[21]=paf_hap1_MHC[[7, 8]].min(axis=1)
paf_hap1_MHC[22]=paf_hap1_MHC[[7, 8]].max(axis=1)

# Sort by contig and start
paf_hap1_MHC = paf_hap1_MHC.sort_values([5, 21], ascending=[True, True])

# Create the "kept intervals"
paf_hap1_MHC_filtered = []

# And now run by row
for idx, row in paf_hap1_MHC.iterrows():
        overlap_found = False           # Default
        for i,k in enumerate(paf_hap1_MHC_filtered):
                if k[5] == row[5] and not (
                        row[22] <= k[21] or 
                        row[21] >= k[22]
                ):
                        # Overlap with a previously kept interval
                        overlap_found = True
                        continue
        if not overlap_found:
                paf_hap1_MHC_filtered.append(row)

paf_hap1_MHC_filtered = pd.DataFrame(paf_hap1_MHC_filtered)

paf_hap1_MHC=paf_hap1_MHC_filtered

paf_hap1_MHC.to_csv(output_paf, sep="\t", index=False, header=None)

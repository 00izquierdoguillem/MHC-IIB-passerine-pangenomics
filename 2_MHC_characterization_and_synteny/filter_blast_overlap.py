import pandas as pd
import re
import sys

if len(sys.argv) != 5:
    print("Usage: python filter_blast_overlap.py <infile> <genes> <outfile> <outfile2>")
    sys.exit(1)

infile = sys.argv[1]
ingenes = sys.argv[2]
outfile = sys.argv[3]
outfile2 = sys.argv[4]

blastable = pd.read_csv(infile, sep="\t", header=None)
genetable = pd.read_csv(ingenes, sep="\t", header=None)

# 1. If gene (col14) in genelist, then add "fastOMA_ortholog" to column15, otherwise write "NA"

blastable[15] = blastable[14].isin(genetable[0]).map(
    {True: "fastOMA_ortholog", False: "NA"}
)

blastable.to_csv(outfile, sep="\t", index=False, header=None)

# 2. Now we want to keep just 1 line per blast hit, prioritizing the ones that have a "fastOMA_ortholog" gene

# Sort by contig,start and ortholog column
blastable = blastable.sort_values([3, 8, 15], ascending=[True, True, False])

blastable_filtered = blastable.drop_duplicates(subset=[3, 8], keep="first")                             # As we have sorted it so that "fastOMA_orthologs" appear first, we will keep only those 1st rows per BLAST hit 

blastable_filtered=pd.DataFrame(blastable_filtered)

blastable_filtered.to_csv(outfile2, sep="\t", index=False, header=None)

print("Reading infile:", infile)
print("Reading genes:", ingenes)
print("Writing outfile1:", outfile)
print("Writing outfile2:", outfile2)

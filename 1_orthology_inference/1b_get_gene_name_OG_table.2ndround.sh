#!/bin/bash
#SBATCH -p edwards,shared # partition to submit to 
#SBATCH -c 1      # cores requested 
#SBATCH -t 0-01:00  # time requested in days-hours:minutes
#SBATCH -o logs/4_get_gene_names_%a.out # name of output file
#SBATCH -e logs/4_get_gene_names_%a.err  # name of error file
#SBATCH --mem=1000 #Mb
#SBATCH --mail-type=END # email me whe nmy job is done running 
#SBATCH --array=1-37%10

WORKDIR=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/2_MHC_project/fastOMA

cd $WORKDIR

KELSIE=/n/netscratch/edwards_lab/Lab/kelsielopez/MHC_test_scott/fastoma
FASTOMA=$KELSIE/out_folder_5_species_conda

GENE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" gene.list)                                      # To run for every previously-ID'd gene

#GENE="RXRB"

# 1. First search gene in isoforms FOR CHICKEN

grep "^$GENE"$'\t' $KELSIE/isoforms_galGal.tsv | cut -f 2 | sed 's/rna-//g' | cut -f 1 -d "." > $GENE.isoforms.tsv

while IFS= read -r read; do

        # 2. Then search the preaviously ID'd isoforms in the fastOMA results and get the HOGs

        cat $FASTOMA/RootHOGs.tsv | grep $read | cut -f 1 > $read.HOGs.tsv
        while IFS= read -r HOG; do
                
                # 3. Search all reads that belongs to this HOG (and create read file for each HOG file if it hasn't been created before)

                if [ ! -f "$HOG.isoforms.tsv" ]; then
                cat $FASTOMA/RootHOGs.tsv | grep $HOG | cut -f 2 | sed 's/rna-//g' | cut -f 1 -d "." > $HOG.isoforms.tsv
                while IFS= read -r final_read; do

                        # 4. And finally get the gene-names for said reads in the isoform files of all species
                        grep $final_read $KELSIE/isoforms_HouseFinch.tsv | awk -v hog="$HOG" -v species="haeMex" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.haeMex.tsv
                        grep $final_read $KELSIE/isoforms_galGal.tsv | awk -v hog="$HOG" -v species="galGal" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.galGal.tsv
                        grep $final_read $KELSIE/isoforms_taeGut.tsv | awk -v hog="$HOG" -v species="taeGut" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.taeGut.tsv
                        grep $final_read $KELSIE/isoforms_chiLan.tsv | awk -v hog="$HOG" -v species="chiLan" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.chiLan.tsv
                        grep $final_read $KELSIE/isoforms_serCan.tsv | awk -v hog="$HOG" -v species="serCan" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.serCan.tsv
                        grep $final_read $KELSIE/isoforms_HouseFinch.tsv | awk -v hog="$HOG" -v species="haeMex" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.haeMex.tsv
                done < $HOG.isoforms.tsv
                fi

        done < $read.HOGs.tsv

done < $GENE.isoforms.tsv


# Repeat for the 3 other species (previous was genes identified in chicken)

grep "^$GENE"$'\t' $KELSIE/isoforms_taeGut.tsv | cut -f 2 | sed 's/rna-//g' | cut -f 1 -d "." > $GENE.isoforms.tsv

while IFS= read -r read; do

        # 2. Then search the preaviously ID'd isoforms in the fastOMA results and get the HOGs

        cat $FASTOMA/RootHOGs.tsv | grep $read | cut -f 1 > $read.HOGs.tsv
        while IFS= read -r HOG; do

                # 3. Search all reads that belongs to this HOG (and create read file for each HOG file if it hasn't been created before)

                if [ ! -f "$HOG.isoforms.tsv" ]; then
                cat $FASTOMA/RootHOGs.tsv | grep $HOG | cut -f 2 | sed 's/rna-//g' | cut -f 1 -d "." > $HOG.isoforms.tsv
                while IFS= read -r final_read; do

                        # 4. And finally get the gene-names for said reads in the isoform files of all species

                        grep $final_read $KELSIE/isoforms_galGal.tsv | awk -v hog="$HOG" -v species="galGal" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.galGal.tsv
                        grep $final_read $KELSIE/isoforms_taeGut.tsv | awk -v hog="$HOG" -v species="taeGut" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.taeGut.tsv
                        grep $final_read $KELSIE/isoforms_chiLan.tsv | awk -v hog="$HOG" -v species="chiLan" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.chiLan.tsv
                        grep $final_read $KELSIE/isoforms_serCan.tsv | awk -v hog="$HOG" -v species="serCan" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.serCan.tsv
                        grep $final_read $KELSIE/isoforms_HouseFinch.tsv | awk -v hog="$HOG" -v species="haeMex" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.haeMex.tsv
                done < $HOG.isoforms.tsv
                fi

        done < $read.HOGs.tsv

done < $GENE.isoforms.tsv

#---

grep "^$GENE"$'\t' $KELSIE/isoforms_serCan.tsv | cut -f 2 | sed 's/rna-//g' | cut -f 1 -d "." > $GENE.isoforms.tsv

while IFS= read -r read; do

        # 2. Then search the preaviously ID'd isoforms in the fastOMA results and get the HOGs

        cat $FASTOMA/RootHOGs.tsv | grep $read | cut -f 1 > $read.HOGs.tsv
        while IFS= read -r HOG; do

                # 3. Search all reads that belongs to this HOG (and create read file for each HOG file if it hasn't been created before)

                if [ ! -f "$HOG.isoforms.tsv" ]; then
                cat $FASTOMA/RootHOGs.tsv | grep $HOG | cut -f 2 | sed 's/rna-//g' | cut -f 1 -d "." > $HOG.isoforms.tsv
                while IFS= read -r final_read; do

                        # 4. And finally get the gene-names for said reads in the isoform files of all species

                        grep $final_read $KELSIE/isoforms_galGal.tsv | awk -v hog="$HOG" -v species="galGal" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.galGal.tsv
                        grep $final_read $KELSIE/isoforms_taeGut.tsv | awk -v hog="$HOG" -v species="taeGut" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.taeGut.tsv
                        grep $final_read $KELSIE/isoforms_chiLan.tsv | awk -v hog="$HOG" -v species="chiLan" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.chiLan.tsv
                        grep $final_read $KELSIE/isoforms_serCan.tsv | awk -v hog="$HOG" -v species="serCan" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.serCan.tsv
                        grep $final_read $KELSIE/isoforms_HouseFinch.tsv | awk -v hog="$HOG" -v species="haeMex" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.haeMex.tsv
                done < $HOG.isoforms.tsv
                fi

        done < $read.HOGs.tsv

done < $GENE.isoforms.tsv


# ---

grep "^$GENE"$'\t' $KELSIE/isoforms_chiLan.tsv | cut -f 2 | sed 's/rna-//g' | cut -f 1 -d "." > $GENE.isoforms.tsv

while IFS= read -r read; do

        # 2. Then search the preaviously ID'd isoforms in the fastOMA results and get the HOGs

        cat $FASTOMA/RootHOGs.tsv | grep $read | cut -f 1 > $read.HOGs.tsv
        while IFS= read -r HOG; do

                # 3. Search all reads that belongs to this HOG (and create read file for each HOG file if it hasn't been created before)

                if [ ! -f "$HOG.isoforms.tsv" ]; then
                cat $FASTOMA/RootHOGs.tsv | grep $HOG | cut -f 2 | sed 's/rna-//g' | cut -f 1 -d "." > $HOG.isoforms.tsv
                while IFS= read -r final_read; do

                        # 4. And finally get the gene-names for said reads in the isoform files of all species

                        grep $final_read $KELSIE/isoforms_galGal.tsv | awk -v hog="$HOG" -v species="galGal" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.galGal.tsv
                        grep $final_read $KELSIE/isoforms_taeGut.tsv | awk -v hog="$HOG" -v species="taeGut" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.taeGut.tsv
                        grep $final_read $KELSIE/isoforms_chiLan.tsv | awk -v hog="$HOG" -v species="chiLan" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.chiLan.tsv
                        grep $final_read $KELSIE/isoforms_serCan.tsv | awk -v hog="$HOG" -v species="serCan" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.serCan.tsv
                        grep $final_read $KELSIE/isoforms_HouseFinch.tsv | awk -v hog="$HOG" -v species="haeMex" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.haeMex.tsv
                done < $HOG.isoforms.tsv
                fi

        done < $read.HOGs.tsv

done < $GENE.isoforms.tsv

# ----

grep "^$GENE"$'\t' $KELSIE/isoforms_HouseFinch.tsv | cut -f 2 | sed 's/rna-//g' | cut -f 1 -d "." > $GENE.isoforms.tsv

while IFS= read -r read; do

        # 2. Then search the preaviously ID'd isoforms in the fastOMA results and get the HOGs

        cat $FASTOMA/RootHOGs.tsv | grep $read | cut -f 1 > $read.HOGs.tsv
        while IFS= read -r HOG; do

                # 3. Search all reads that belongs to this HOG (and create read file for each HOG file if it hasn't been created before)

                if [ ! -f "$HOG.isoforms.tsv" ]; then
                cat $FASTOMA/RootHOGs.tsv | grep $HOG | cut -f 2 | sed 's/rna-//g' | cut -f 1 -d "." > $HOG.isoforms.tsv
                while IFS= read -r final_read; do

                        # 4. And finally get the gene-names for said reads in the isoform files of all species

                        grep $final_read $KELSIE/isoforms_galGal.tsv | awk -v hog="$HOG" -v species="galGal" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.galGal.tsv
                        grep $final_read $KELSIE/isoforms_taeGut.tsv | awk -v hog="$HOG" -v species="taeGut" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.taeGut.tsv
                        grep $final_read $KELSIE/isoforms_chiLan.tsv | awk -v hog="$HOG" -v species="chiLan" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.chiLan.tsv
                        grep $final_read $KELSIE/isoforms_serCan.tsv | awk -v hog="$HOG" -v species="serCan" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.serCan.tsv
                        grep $final_read $KELSIE/isoforms_HouseFinch.tsv | awk -v hog="$HOG" -v species="haeMex" 'BEGIN{OFS="\t"} {print $0, hog, species}' >> $HOG.gene_names.haeMex.tsv
                done < $HOG.isoforms.tsv
                fi

        done < $read.HOGs.tsv

done < $GENE.isoforms.tsv

#############

# Finally join all gene_names into a single file

cat *gene_names.*.tsv > fastOMA.MHC_rootHOGs.tsv
sort fastOMA.MHC_rootHOGs.tsv | uniq > prova
mv prova fastOMA.MHC_rootHOGs.tsv

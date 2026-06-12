library(tidyr)
library(ape)
library(pegas)
library(PopGenome)
library(ggplot2)
library(ggpattern)
library(grid)
library(patchwork)

summary_table_global<- read.table("allspecies.final_blast_summary.exon2_exon3_combined.Feb26.Scott.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE,comment.char = "")
summary_table_global <- summary_table_global[summary_table_global$SPECIES!="CY" & summary_table_global$SPECIES!="AA",]
summary_table_global[summary_table_global$SPECIES=="HF_E",]$SPECIES <- "HF"
summary_table_global[summary_table_global$SPECIES=="HF_W",]$SPECIES <- "HF"
summary_table_global[summary_table_global$CONTIG_TYPE=="Uncertain",]$CONTIG_TYPE <- "Core" 			# Do this if you consider uncertain contigs as core

alignments <- list(
  ex2_SJ = read.dna("allspecies_outgroups.ex2_mafft.SJ_def.fas", format="fasta"),
  ex2_HF = read.dna("allspecies_outgroups.ex2_mafft.HF_def.fas", format="fasta"),
  ex2_HM = read.dna("allspecies_outgroups.ex2_mafft.HM_def.fas", format="fasta"),
  ex3_SJ = read.dna("allspecies_outgroups_ex3_mafft.SJ_def.fas", format="fasta"),
  ex3_HF = read.dna("allspecies_outgroups_ex3_mafft.HF_def.fas", format="fasta"),
  ex3_HM = read.dna("allspecies_outgroups_ex3_mafft.HM_def.fas", format="fasta")
)

compute_pi_table <- function(alignments, summary_table) {
  
  species_vec <- unique(summary_table$SPECIES)
  contig_vec  <- unique(summary_table$GENE_TYPE)					# Create species and contig_type lists
  
  results_list <- list()											# Initialize list and counter
  counter <- 1
  
  for (exon_name in names(alignments)) {							# You consider the alignment name to be the exon name (?) - iterate through files
    
    dna <- alignments[[exon_name]]
    
    for (sp in species_vec) {											# Now loop across species and contig_types
      for (ct in contig_vec) {
        
        df_subset <- summary_table[										# Subset the metadata table with species & contig type information
          summary_table$SPECIES == sp &
          summary_table$GENE_TYPE == ct,
        ]
        
        dna_subset <- dna[rownames(dna) %in% df_subset$BLAST_ID, ]		# Now subset the fasta alignment to keep only BLAST hits that coincide with the metadata filter
        
        if (nrow(dna_subset) < 2) {
          pi_value <- NA
        } else {
          pi_value <- nuc.div(dna_subset)								# Calculate pi
        }
        
        results_list[[counter]] <- data.frame(							# Output the pi (and species, contig and exon information) to a row for a dataframe. The counter is the row index, and it will keep adding according to the next line of code (creating a list of rows)
          SPECIES = sp,
          GENE_TYPE = ct,
          EXON = exon_name,
          PI = pi_value
        )
        
        counter <- counter + 1
      }
    }
  }
  
  do.call(rbind, results_list)											# Transform to df
}

## Run function

pi_table <- compute_pi_table(alignments, summary_table_global)

pi_table[pi_table$EXON=="ex2_SJ" | pi_table$EXON=="ex2_HF" | pi_table$EXON=="ex2_HM" ,]$EXON <- "exon2"
pi_table[pi_table$EXON=="ex3_SJ" | pi_table$EXON=="ex3_HF" | pi_table$EXON=="ex3_HM" ,]$EXON <- "exon3"				# Rename exon column

#################### Now plot

palette_SVE <- c(
  AA   = "#F7F00B",
  AC   = "#A9ADB8",
  AI   = "#114EF2",
  AW   = "#8AA5EE",
  HM = "#8BB19C",
  HF_E  = "#EF3B39",
  HF_W = "#EF9D39",
  HF  = "#EF3B39",
  CY   = "#231F20")

pi_table$SPECIES <- factor(pi_table$SPECIES,  levels=c("AI","AW","AC","HF","HM"),ordered = TRUE)

a <- ggplot(pi_table[pi_table$EXON=="exon2",], aes(x = SPECIES, y = PI, fill = SPECIES, pattern=GENE_TYPE)) +
  geom_col_pattern(aes(pattern = GENE_TYPE, pattern_fill = SPECIES), color = "black",pattern_fill = "black",pattern_angle = 45,pattern_spacing = 0.02,pattern_density = 0.05, position="dodge") +
   labs(y = "Nucleotide diversity in exon2 (bp⁻¹)",
    x = "Species") +
  theme_bw() + scale_pattern_manual(values = c(Potential_pseudogene = "crosshatch",Potential_coding_gene = "none")) +
  scale_fill_manual(values=palette_SVE) + scale_y_continuous(limits=c(-0.002, 0.18))

a


b <- ggplot(pi_table[pi_table$EXON=="exon3",], aes(x = SPECIES, y = PI, fill = SPECIES, pattern=GENE_TYPE)) +
  geom_col_pattern(aes(pattern = GENE_TYPE, pattern_fill = SPECIES), color = "black",pattern_fill = "black",pattern_angle = 45,pattern_spacing = 0.02,pattern_density = 0.05, position="dodge") +
   labs(y = "Nucleotide diversity in exon3 (bp⁻¹)",
    x = "Species") +
  theme_bw() + scale_pattern_manual(values = c(Potential_pseudogene = "crosshatch",Potential_coding_gene = "none")) +
  scale_fill_manual(values=palette_SVE) + scale_y_continuous(limits=c(-0.002, 0.18))

b

combined_plot <- a + b + plot_layout(ncol = 2, guides= "collect")
combined_plot
ggsave ("../nuc_div/MHC.nuc_div.pdf", combined_plot)

#########################################################################################################################
#################### Now do the same to get nucleotide diversity per individual #####################################################################
#########################################################################################################################

compute_pi_ind <- function(alignments, summary_table) {
  
  species_vec <- unique(summary_table$SPECIES)
  contig_vec  <- unique(summary_table$GENE_TYPE)					# Create species and contig_type lists
  ind_vec <- unique(summary_table$IND)

  results_list <- list()											# Initialize list and counter
  counter <- 1
  
  for (exon_name in names(alignments)) {							# You consider the alignment name to be the exon name (?) - iterate through files
    
    dna <- alignments[[exon_name]]
    
    for (sp in species_vec) {											# Now loop across species and contig_types
      for (ct in contig_vec) {
      for(ind in ind_vec) {  
        df_subset <- summary_table[										# Subset the metadata table with species & contig type information
          summary_table$SPECIES == sp &
          summary_table$GENE_TYPE == ct &
          summary_table$IND == ind,
        ]
        
        dna_subset <- dna[rownames(dna) %in% df_subset$BLAST_ID, ]		# Now subset the fasta alignment to keep only BLAST hits that coincide with the metadata filter
        
        if (nrow(dna_subset) < 2) {
          pi_value <- NA
        } else {
          pi_value <- nuc.div(dna_subset)								# Calculate pi
        }
        
        results_list[[counter]] <- data.frame(							# Output the pi (and species, contig and exon information) to a row for a dataframe. The counter is the row index, and it will keep adding according to the next line of code (creating a list of rows)
          SPECIES = sp,
          GENE_TYPE = ct,
          EXON = exon_name,
          IND= ind,
          PI = pi_value
        )
        
        counter <- counter + 1
      }
  	}
    }
  }
  
  do.call(rbind, results_list)											# Transform to df
}

## Run function

pi_table_ind <- compute_pi_ind(alignments, summary_table_global)


pi_table_ind[pi_table_ind$EXON=="ex2_SJ" | pi_table_ind$EXON=="ex2_HF" | pi_table_ind$EXON=="ex2_HM" ,]$EXON <- "exon2"
pi_table_ind[pi_table_ind$EXON=="ex3_SJ" | pi_table_ind$EXON=="ex3_HF" | pi_table_ind$EXON=="ex3_HM" ,]$EXON <- "exon3"				# Rename exon column

#################### Now plot

pi_table_ind$SPECIES <- factor(pi_table_ind$SPECIES,  levels=c("AI","AW","AC","HF","HM"),ordered = TRUE)

c <- ggplot(pi_table_ind[pi_table_ind$EXON=="exon2",], aes(x = SPECIES, y = PI, fill = SPECIES, pattern=GENE_TYPE)) +
  geom_boxplot_pattern(aes(pattern = GENE_TYPE, pattern_fill = SPECIES), color = "black",pattern_fill = "black",pattern_angle = 45,pattern_spacing = 0.02,pattern_density = 0.05) +
   labs(y = "Nucleotide diversity in exon2 per individual (bp⁻¹)",
    x = "Species") +
  theme_bw() + scale_pattern_manual(values = c(Potential_pseudogene = "crosshatch",Potential_coding_gene = "none")) +
  scale_fill_manual(values=palette_SVE) + scale_y_continuous(limits=c(-0.002, 0.28))

c


d <- ggplot(pi_table_ind[pi_table_ind$EXON=="exon3",], aes(x = SPECIES, y = PI, fill = SPECIES, pattern=GENE_TYPE)) +
  geom_boxplot_pattern(aes(pattern = GENE_TYPE, pattern_fill = SPECIES), color = "black",pattern_fill = "black",pattern_angle = 45,pattern_spacing = 0.02,pattern_density = 0.05) +
   labs(y = "Nucleotide diversity in exon3 per individual (bp⁻¹)",
    x = "Species") +
  theme_bw() + scale_pattern_manual(values = c(Potential_pseudogene = "crosshatch",Potential_coding_gene = "none")) +
  scale_fill_manual(values=palette_SVE) + scale_y_continuous(limits=c(-0.002, 0.28))

d

combined_plot <- c + d + plot_layout(ncol = 2, guides= "collect")
combined_plot
ggsave ("../nuc_div/MHC.nuc_div.ind.pdf", combined_plot)



#########################################################################################################################
#################### Now do the same to get diversity per haplotype #####################################################################
#########################################################################################################################

compute_pi_hap <- function(alignments, summary_table) {
  
  species_vec <- unique(summary_table$SPECIES)
  contig_vec  <- unique(summary_table$GENE_TYPE)					# Create species and contig_type lists
  hap_vec <- unique(summary_table$HAP)

  results_list <- list()											# Initialize list and counter
  counter <- 1
  
  for (exon_name in names(alignments)) {							# You consider the alignment name to be the exon name (?) - iterate through files
    
    dna <- alignments[[exon_name]]
    
    for (sp in species_vec) {											# Now loop across species and contig_types
      for (ct in contig_vec) {
      for(hap in hap_vec) {  
        df_subset <- summary_table[										# Subset the metadata table with species & contig type information
          summary_table$SPECIES == sp &
          summary_table$GENE_TYPE == ct &
          summary_table$HAP == hap,
        ]
        
        dna_subset <- dna[rownames(dna) %in% df_subset$BLAST_ID, ]		# Now subset the fasta alignment to keep only BLAST hits that coincide with the metadata filter
        
        if (nrow(dna_subset) < 2) {
          pi_value <- NA
        } else {
          pi_value <- nuc.div(dna_subset)								# Calculate pi
        }
        
        results_list[[counter]] <- data.frame(							# Output the pi (and species, contig and exon information) to a row for a dataframe. The counter is the row index, and it will keep adding according to the next line of code (creating a list of rows)
          SPECIES = sp,
          GENE_TYPE = ct,
          EXON = exon_name,
          HAP= hap,
          PI = pi_value
        )
        
        counter <- counter + 1
      }
  	}
    }
  }
  
  do.call(rbind, results_list)											# Transform to df
}

## Run function

pi_table_hap <- compute_pi_hap(alignments, summary_table_global)


pi_table_hap[pi_table_hap$EXON=="ex2_SJ" | pi_table_hap$EXON=="ex2_HF" | pi_table_hap$EXON=="ex2_HM" ,]$EXON <- "exon2"
pi_table_hap[pi_table_hap$EXON=="ex3_SJ" | pi_table_hap$EXON=="ex3_HF" | pi_table_hap$EXON=="ex3_HM" ,]$EXON <- "exon3"				# Rename exon column

#################### Now plot

pi_table_hap$SPECIES <- factor(pi_table_hap$SPECIES,  levels=c("AI","AW","AC","HF","HM"),ordered = TRUE)

e <- ggplot(pi_table_hap[pi_table_hap$EXON=="exon2",], aes(x = SPECIES, y = PI, fill = SPECIES, pattern=GENE_TYPE)) +
  geom_boxplot_pattern(aes(pattern = GENE_TYPE, pattern_fill = SPECIES), color = "black",pattern_fill = "black",pattern_angle = 45,pattern_spacing = 0.02,pattern_density = 0.05) +
   labs(y = "Nucleotide diversity in exon2 per haplotype",
    x = "Species") +
  theme_bw() + scale_pattern_manual(values = c(Potential_pseudogene = "crosshatch",Potential_coding_gene = "none")) +
  scale_fill_manual(values=palette_SVE) + scale_y_continuous(limits=c(-0.002, 0.27))

e


f <- ggplot(pi_table_hap[pi_table_hap$EXON=="exon3",], aes(x = SPECIES, y = PI, fill = SPECIES, pattern=GENE_TYPE)) +
  geom_boxplot_pattern(aes(pattern = GENE_TYPE, pattern_fill = SPECIES), color = "black",pattern_fill = "black",pattern_angle = 45,pattern_spacing = 0.02,pattern_density = 0.05) +
   labs(y = "Nucleotide diversity in exon3 per haplotype",
    x = "Species") +
  theme_bw() + scale_pattern_manual(values = c(Potential_pseudogene = "crosshatch",Potential_coding_gene = "none")) +
  scale_fill_manual(values=palette_SVE) + scale_y_continuous(limits=c(-0.002, 0.27), transform="sqrt")

f

combined_plot <- e + f + plot_layout(ncol = 2, guides= "collect")
combined_plot
ggsave ("../nuc_div/MHC.nuc_div.hap.pdf", combined_plot)

#------------------------

combined_plot <- a+b+c+d + plot_layout(ncol = 2, guides= "collect")
combined_plot
ggsave ("../nuc_div/MHC.nuc_div.allHF.pdf", combined_plot)

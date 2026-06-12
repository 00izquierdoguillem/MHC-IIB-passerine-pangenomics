library(tidyr)
library(ape)
library(ggplot2)
library(ggpattern)
library(grid)
library(patchwork)

## This script should allow you to get nucleotid diversity density plots

summary_table_global<- read.table("allspecies.final_blast_summary.exon2_exon3_combined.Feb26.Scott.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE,comment.char = "")
summary_table_global <- summary_table_global[summary_table_global$SPECIES!="CY" & summary_table_global$SPECIES!="AA",]
summary_table_global[summary_table_global$SPECIES=="HF_E",]$SPECIES <- "HF"
summary_table_global[summary_table_global$SPECIES=="HF_W",]$SPECIES <- "HF"
#summary_table_global[summary_table_global$CONTIG_TYPE=="Uncertain",]$CONTIG_TYPE <- "Core"       # Do this if you consider uncertain contigs as core

alignments <- list(
  ex2_SJ = read.dna("allspecies_outgroups.ex2_mafft.SJ_def.fas", format="fasta"),
  ex2_HF = read.dna("allspecies_outgroups.ex2_mafft.HF_def.fas", format="fasta"),
  ex2_HM = read.dna("allspecies_outgroups.ex2_mafft.HM_def.fas", format="fasta"),
  ex3_SJ = read.dna("allspecies_outgroups_ex3_mafft.SJ_def.fas", format="fasta"),
  ex3_HF = read.dna("allspecies_outgroups_ex3_mafft.HF_def.fas", format="fasta"),
  ex3_HM = read.dna("allspecies_outgroups_ex3_mafft.HM_def.fas", format="fasta")
)

# 1. Function that differentiates comparisons and returns distances

split_dist_matrix_vec <- function(dna_alignment, summary_table) {
  
  # Distance matrix
  dist_mat <- as.matrix(dist.dna(dna_alignment))
  
  ## Remove self-comparisons
  diag(dist_mat) <- NA

  # Sequence names
  seq_names <- rownames(dist_mat)
  
  # Map contigs to contig type
  contig_type <- summary_table$CONTIG_TYPE[match(seq_names, summary_table$BLAST_ID)]
  names(contig_type) <- seq_names
  
  # Logical matrices describing pair categories (that is, masking the values that don't fit the condition)
  core_i <- contig_type == "Core"
  notcore_i <- contig_type == "Not_core"
  
  core_core_mask <- outer(core_i, core_i, "&")
  notcore_notcore_mask <- outer(notcore_i, notcore_i, "&")
  core_notcore_mask <- outer(core_i, notcore_i, "&") | outer(notcore_i, core_i, "&")
  
  # Create matrices
  core_core <- dist_mat
  notcore_notcore <- dist_mat
  core_notcore <- dist_mat
  
  core_core[!core_core_mask] <- NA
  notcore_notcore[!notcore_notcore_mask] <- NA
  core_notcore[!core_notcore_mask] <- NA
  
  return(list(
    core_core = core_core,
    notcore_notcore = notcore_notcore,
    core_notcore = core_notcore
  ))
}

# 2. Apply the function

dist_split <- lapply(alignments, split_dist_matrix_vec, 
                     summary_table = summary_table_global)

# 3. Transform to a dataframe with columns describing exon, species and comparison type.

build_distance_dataframe_vec <- function(alignments, summary_table){
  
  dist_split <- lapply(alignments, split_dist_matrix_vec,
                       summary_table = summary_table)
  
  all_dfs <- lapply(names(dist_split), function(aln_name){
    
    mats <- dist_split[[aln_name]]
    
    exon <- ifelse(grepl("ex2", aln_name), "exon2", "exon3")
    species <- sub(".*_", "", aln_name)
    
    dfs <- lapply(names(mats), function(comp_type){
      
      mat <- mats[[comp_type]]
      
      df <- as.data.frame(as.table(mat))
      colnames(df) <- c("SEQ1","SEQ2","DISTANCE")
      
      df$SPECIES <- species
      df$EXON <- exon
      df$COMPARISON_TYPE <- comp_type
      
      df
      
    })
    
    do.call(rbind, dfs)
    
  })
  
  dist_matrix_global <- do.call(rbind, all_dfs)
  
  return(dist_matrix_global)
}

dist_matrix_global <- build_distance_dataframe_vec(
  alignments,
  summary_table_global
)

# 4. Remove NAs

dist_matrix_global <- dist_matrix_global[!is.na(dist_matrix_global$DISTANCE),]

# 5. Plot

c <- ggplot(dist_matrix_global[dist_matrix_global$EXON=="exon2" & dist_matrix_global$COMPARISON_TYPE=="core_core",], aes(x = DISTANCE, fill = SPECIES, colour=SPECIES)) +
  geom_density(alpha=0.3, adjust=1.5) +
  theme_classic() + theme(legend.position="none") +scale_fill_manual(values=palette_SVE) +
  labs(y = "Density",
    x = "Genetic distance between \n core MHC-IIB exon 2 copies") + scale_x_continuous(limits=c(0,0.5))

c

d <- ggplot(dist_matrix_global[dist_matrix_global$EXON=="exon2" & dist_matrix_global$COMPARISON_TYPE=="notcore_notcore",], aes(x = DISTANCE, fill = SPECIES, colour=SPECIES)) +
  geom_density_pattern(aes(pattern = COMPARISON_TYPE, pattern_fill = SPECIES), pattern_fill = "black",pattern_angle = 45,pattern_spacing = 0.02,pattern_density = 0.05,
    alpha=0.3, adjust=2.5, pattern_alpha=0.05) +
  scale_fill_manual(values=palette_SVE) + scale_pattern_manual(values = c("notcore_notcore" = "crosshatch")) + theme_classic() + theme(legend.position="none") +
  labs(y = "Density",
    x = "Genetic distance between \n outlying MHC-IIB exon 2 copies") + scale_x_continuous(limits=c(0,0.5))

d

e <- ggplot(dist_matrix_global[dist_matrix_global$EXON=="exon2" & dist_matrix_global$COMPARISON_TYPE=="core_notcore",], aes(x = DISTANCE, fill = SPECIES, colour=SPECIES)) +
  geom_density_pattern(aes(pattern = COMPARISON_TYPE, pattern_fill = SPECIES), pattern_fill = "black",pattern_angle = 45,pattern_spacing = 0.02,pattern_density = 0.05,
    alpha=0.3, adjust=2, pattern_alpha=0.05) +
  scale_fill_manual(values=palette_SVE) + scale_pattern_manual(values = c("core_notcore" = "stripe")) + theme_classic() + theme(legend.position="none") +
  labs(y = "Density",
    x = "Genetic distance between core and \n outlying MHC-IIB exon 2 copies") + scale_x_continuous(limits=c(0,0.5))

e

f <- ggplot(dist_matrix_global[dist_matrix_global$EXON=="exon3" & dist_matrix_global$COMPARISON_TYPE=="core_core",], aes(x = DISTANCE, fill = SPECIES, colour=SPECIES)) +
  geom_density(alpha=0.3, adjust=3) +
  theme_bw() +  scale_fill_manual(values=palette_SVE)

f

g <- ggplot(dist_matrix_global[dist_matrix_global$EXON=="exon3" & dist_matrix_global$COMPARISON_TYPE=="notcore_notcore",], aes(x = DISTANCE, fill = SPECIES, colour=SPECIES)) +
  geom_density(alpha=0.3, adjust=4) +
  theme_bw() +  scale_fill_manual(values=palette_SVE)

g

h <- ggplot(dist_matrix_global[dist_matrix_global$EXON=="exon3" & dist_matrix_global$COMPARISON_TYPE=="core_notcore",], aes(x = DISTANCE, fill = SPECIES, colour=SPECIES)) +
  geom_density(alpha=0.3, adjust=3) +
  theme_bw() +  scale_fill_manual(values=palette_SVE)

h

###### Prova tot junt

k <- ggplot(dist_matrix_global[dist_matrix_global$EXON=="exon2" ,], aes(x = DISTANCE, fill = SPECIES, colour=SPECIES, pattern=COMPARISON_TYPE)) +
  geom_density_pattern(aes(pattern = COMPARISON_TYPE, pattern_fill = SPECIES), pattern_fill = "black",pattern_angle = 45,pattern_spacing = 0.02,pattern_density = 0.05,
    alpha=0.3, adjust=2.5, pattern_alpha=0.05) +
  theme_bw() +  scale_fill_manual(values=palette_SVE) + scale_pattern_manual(values = c("notcore_notcore" = "crosshatch", "core_core" = "none", "core_notcore" = "stripe"))

k

library(ape)
library(phangorn)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# 1. Exon2 dist matrix

exon2_tree <- read.tree(file = "allspecies_outgroups.ex2_mafft.fa.tre")
exon2_dist <- cophenetic(exon2_tree)
exon2_dist <- as.data.frame(exon2_dist)

exon2_names <- exon2_tree$tip.label
exon2_names <- as.data.frame(exon2_names)

rownames(exon2_dist) <- exon2_names$exon2_names

### 2. Now exon3

exon3_tree <- read.tree(file = "allspecies_outgroups_ex3_mafft.fa.tre")
exon3_dist <- cophenetic(exon3_tree)
exon3_dist <- as.data.frame(exon3_dist)

exon3_names <- exon3_tree$tip.label
exon3_names <- as.data.frame(exon3_names)

rownames(exon3_dist) <- exon3_names$exon3_names

## 3. And finally load the matrix from the "BLAST results with frameshifts" matrix to create the distance matrix

blast_results <- read.table("allspecies.final_blast_summary.exon2_exon3_combined.Feb26.csv", sep="\t", header=TRUE, comment.char="")
blast_results[blast_results$CONTIG_TYPE=="Uncertain",]$CONTIG_TYPE <- "Core" 
start_diff <- abs(outer(blast_results$START, blast_results$START, "-"))

# Logical matrices for matching CONTIG and EXON
same_contig <- outer(blast_results$CONTIG, blast_results$CONTIG, "==")
same_exon   <- outer(blast_results$EXON, blast_results$EXON, "==")

# Keep only valid comparisons
dist_matrix <- start_diff
dist_matrix[!(same_contig & same_exon)] <- NA

# Add row/column names
rownames(dist_matrix) <- blast_results$BLAST_ID
colnames(dist_matrix) <- blast_results$BLAST_ID

dist_matrix<- as.data.frame(dist_matrix)

# 4. Compare dist matrices

# Create long versions of each matrix and then create a new column with pairwise comparison name. Drop NAs

dist_matrix_long <- dist_matrix %>% rownames_to_column("Hit1") %>%
  pivot_longer(
    cols = -Hit1,
    names_to = "Hit2",
    values_to = "Phys_dist"
  )

dist_matrix_long <- dist_matrix_long[dist_matrix_long$Phys_dist!=0,]					# Remvoe same-hit comparisons or the ones flagged as NA
dist_matrix_long <- dist_matrix_long[!is.na(dist_matrix_long$Phys_dist),]

dist_matrix_long$comparison_name<- paste0(dist_matrix_long$Hit1, dist_matrix_long$Hit2)

## Same with phylo dist matrices

exon2_long <- exon2_dist %>% rownames_to_column("Hit1") %>%
  pivot_longer(
    cols = -Hit1,
    names_to = "Hit2",
    values_to = "Phylo_dist"
  )

exon2_long$comparison_name<- paste0(exon2_long$Hit1, exon2_long$Hit2)

exon2_long <- exon2_long[exon2_long$comparison_name %in% dist_matrix_long$comparison_name,]				#### Keep only the previously filtered comparisons

## Exon3 

exon3_long <- exon3_dist %>% rownames_to_column("Hit1") %>%
  pivot_longer(
    cols = -Hit1,
    names_to = "Hit2",
    values_to = "Phylo_dist"
  )

exon3_long$comparison_name<- paste0(exon3_long$Hit1, exon3_long$Hit2)

exon3_long <- exon3_long[exon3_long$comparison_name %in% dist_matrix_long$comparison_name,]

# 5. Join matrices and plot

plot_exon2 <- merge(dist_matrix_long, exon2_long)
plot_exon3 <- merge(dist_matrix_long, exon3_long)

write.table(plot_exon2, "checkpoint_exon2.csv")
write.table(plot_exon3, "checkpoint_exon3.csv")


# Plot

palette_SVE <- c(
  AA   = "#F7F00B",
  AC   = "#A9ADB8",
  AI   = "#114EF2",
  AW   = "#8AA5EE",
  HM = "#8BB19C",
  HF_E  = "#EF3B39",
  HF_W = "#EF9D39",
  HF  = "#EF3B39",
  CY   = "#231F20"
)

a <- ggplot(plot_exon2, aes(x = Phys_dist, y = Phylo_dist)) +
  geom_point(size=0.3, alpha=0.2) +
   labs(y = "Phylogenetic distance",
    x = "Physical distance") +
  theme_bw() #+
#  scale_color_manual(values=palette_SVE)

a

### Option! Colour by contig (core vs non-core!). But maybe better color by species? And 2 separate plots?
### To do so, we have to include metadata here

blast_results <- unique(subset(blast_results, select=c(HAP,SPECIES,GENE_TYPE,CONTIG_TYPE,BLAST_ID)))

blast_results <- blast_results[blast_results$SPECIES!="CY" & blast_results$SPECIES!="AA",]
blast_results[blast_results$SPECIES=="HF_E",]$SPECIES <- "HF"
blast_results[blast_results$SPECIES=="HF_W",]$SPECIES <- "HF"

plot_exon2 <- merge(plot_exon2, blast_results, by.x="Hit1", by.y="BLAST_ID")
plot_exon3 <- merge(plot_exon3, blast_results, by.x="Hit1", by.y="BLAST_ID")
#
#plot_exon2_final <- plot_exon2[plot_exon2$CONTIG_TYPE == "Not_core",] 
#pval_exon2 <- plot_exon2_final[plot_exon2_final$SPECIES=="HF",]
#model_HF <- lm(Phylo_dist ~ Phys_dist, data = pval_exon2)
#summary(model_HF)


b <- ggplot(plot_exon2, aes(x = Phys_dist, y = Phylo_dist, col=SPECIES)) +
  geom_point(size=0.3, alpha=0.2) +
   labs(y = "Phylogenetic distance",
    x = "Physical distance") +
  theme_bw() +
  scale_color_manual(values=palette_SVE)

b

c <- ggplot(plot_exon2, aes(x = Phys_dist, y = Phylo_dist, col=CONTIG_TYPE)) +
  geom_point(size=0.3, alpha=0.2) +
   labs(y = "Phylogenetic distance",
    x = "Physical distance") +
  theme_bw() +
  scale_color_manual(values=c("#282EA1", "#D1563D", "#000000"))

c

### Let's separate between core vs not_core

#plot_exon2[plot_exon2$CONTIG_TYPE=="Uncertain",]$CONTIG_TYPE <- "Core" 
plot_exon2_final <- plot_exon2[plot_exon2$CONTIG_TYPE == "Core",]

d <- ggplot(plot_exon2_final, aes(x = Phys_dist, y = Phylo_dist, col=SPECIES)) +
    labs(y = "Phylogenetic distance",
    x = "Physical distance") +
  theme_bw() +
  scale_color_manual(values=palette_SVE) + geom_smooth(method="lm", alpha=0.2) + geom_point(size=0.3, alpha=0.2)

d

plot_exon2_final <- plot_exon2[plot_exon2$CONTIG_TYPE == "Not_core",]

e <- ggplot(plot_exon2_final, aes(x = Phys_dist, y = Phylo_dist, col=SPECIES)) +
    labs(y = "Phylogenetic distance",
    x = "Physical distance") +
  theme_bw() +
  scale_color_manual(values=palette_SVE) + geom_smooth(method="lm", alpha=0.2) + geom_point(size=0.3, alpha=0.2)

  e


###### Repeat with exon 3


f <- ggplot(plot_exon3, aes(x = Phys_dist, y = Phylo_dist, col=SPECIES)) +
  geom_point(size=0.3, alpha=0.2) +
   labs(y = "Phylogenetic distance",
    x = "Physical distance") +
  theme_bw() +
  scale_color_manual(values=palette_SVE)

f

g <- ggplot(plot_exon3, aes(x = Phys_dist, y = Phylo_dist, col=CONTIG_TYPE)) +
  geom_point(size=0.3, alpha=0.2) +
   labs(y = "Phylogenetic distance",
    x = "Physical distance") +
  theme_bw() +
  scale_color_manual(values=c("#282EA1", "#D1563D", "#000000"))

g

### Let's separate between core vs not_core

#plot_exon3[plot_exon3$CONTIG_TYPE=="Uncertain",]$CONTIG_TYPE <- "Core" 
plot_exon3_final <- plot_exon3[plot_exon3$CONTIG_TYPE == "Core",]

h <- ggplot(plot_exon3_final, aes(x = Phys_dist, y = Phylo_dist, col=SPECIES)) +
    labs(y = "Phylogenetic distance",
    x = "Physical distance") +
  theme_bw() +
  scale_color_manual(values=palette_SVE) + geom_smooth(method="lm", alpha=0.2) + geom_point(size=0.3, alpha=0.2)

h

plot_exon3_final <- plot_exon3[plot_exon3$CONTIG_TYPE == "Not_core",]

i <- ggplot(plot_exon3_final, aes(x = Phys_dist, y = Phylo_dist, col=SPECIES)) +
    labs(y = "Phylogenetic distance",
    x = "Physical distance") +
  theme_bw() +
  scale_color_manual(values=palette_SVE) + geom_smooth(method="lm", alpha=0.2) + geom_point(size=0.3, alpha=0.2)

i

######## Plot together


library(patchwork)

combined_plot <- c+ g + d + h + e + i + plot_layout(ncol = 2, guides= "collect")
combined_plot

ggsave ("Phylo_physical_distance.MHC.pdf", combined_plot)

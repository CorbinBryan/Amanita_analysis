library(ggtree)
library(treeio)
library(ggplot2)
library(dplyr)
library(processx)
library(tidyverse)
library(ggtreeExtra)
library(ape)
setwd("/Users/corbinbryan/Desktop/Amanita_analysis/scripts/R")
treefile <- "../../buscos/astral_sp_trees/astral4_anno2.tree"
tree <- read.astral(treefile)

#freq_data <- data.frame(
#  tree@data$node,
#  tree@data$q1,
#  tree@data$q2,
#  tree@data$q3
#  )

#colnames(freq_data) <- c(
#  "node",
#  "q1",
#  "q2",
#  "q3"
#)

rooted_phylo <- root(
  tree@phylo,
  "Volvo1_AssemblyScaffolds_Repeatmasked.fasta",
  resolve.root = TRUE
)

tree_rooted <- tree
tree_rooted@phylo <- rooted_phylo



freq_long <- tree_rooted@data %>%
  select(node, pp1, pp2, pp3) #%>%
  #pivot_longer(
  #  cols = c(q1, q2, q3),
  #  names_to = "topology",
  #  values_to = "value"
  #)

#freq_long$node <- freq_long$node + 1 

exclude <- c(
  "20045_S349_spades_no_contam.fasta",               
  "20031_S348_spades_no_contam.fasta",               
  "GCA_054369615.1_ASM5436961v1_genomic.fna",      
  "30021_S33_L007",                                   
  "GCA_001691765.1_ASM169176v1_genomic.fna",          
  "30000_S38_L007",
  "NagxHevest_S354_spades_no_contam.fasta",
  "NagxHevestB_S355_spades_no_contam.fasta",
  "SRR36137122_contigs.fasta",
  "Gril1_S352_spades_no_contam.fasta",
  "Frag1_S351_spades_no_contam.fasta",
  "Aus332_S350_spades_no_contam.fasta",
  "NzAUS95_S358_spades_no_contam.fasta",
  "Nes1_S356_spades_no_contam.fasta",
  "Skrap3_S361_spades_no_contam.fasta",            
  "11665_S342_spades_no_contam.fasta",              
  "11663_S340_spades_no_contam.fasta",                
  "11662_no_contam.fasta",                            
  "11669_S346_spades_no_contam.fasta",                
  "11668_S345_spades_no_contam.fasta",                
  "11670_S347_spades_no_contam.fasta",                
  "11666_S343_spades_no_contam.fasta",                
  "11667_S344_spades_no_contam.fasta",                
  "11664_S341_spades_no_contam.fasta",                
  "Sogn5_S362_spades_no_contam.fasta",                
  "Ring1_S359_spades_no_contam.fasta",                
  "Kara3_S353_spades_no_contam.fasta",               
  "Wirz_S363_spades_no_contam.fasta"  
)

write_lines(
  subset(tree_rooted@phylo[["tip.label"]], 
         !(tree_rooted@phylo[["tip.label"]] %in% exclude)),
  "tip_names.txt"
)


system("./process_tip_names.sh")

name_key <- read.csv(
  "tips_species.tsv",
  sep = "\t", 
  header = FALSE, 
  col.names = c("tip", "species")
) #%>% subset((!species %in% exclude))


pies <- nodepie(freq_long, cols = c("pp1", "pp2", "pp3"), 
                color = c("pp1" = "black", "pp2" = "white", "pp3" = "red"))
p <- ggtree(tree_rooted) +
  #geom_tiplab() + 
  theme_tree() + 
  geom_cladelabel(node=93, label='italic("A. muscaria")', parse = TRUE, color = "#b98884") + 
  geom_cladelabel(node=113, label='italic("A. chrysoblema")',parse = TRUE, color = "skyblue") + 
  geom_cladelabel(node=115, label='italic("A. flavivolvata")',parse = TRUE, color = "gold") +
  xlim(0, max(tree_rooted@phylo[["edge.length"]]) * 1.1)
  
j <- p %<+% name_key + geom_tiplab(aes(subset=(!label %in% exclude),
                                  label=paste0("italic('", species, "')"), 
                                  parse = TRUE), parse = TRUE) +


# Add pie charts as insets
final_plot <- inset(j, pies, width = 0.02, height = 0.02) 

print(final_plot) 

ggsave("annotated_tree.svg", plot = final_plot)

### ANI Analysis
setwd("/Users/corbinbryan/Desktop/Amanita_analysis/scripts/R")

library(pheatmap)

library(pheatmap)

# ---- 1. Read file ----
lines <- readLines("../../skani_dist_matrix.txt")
n <- as.integer(lines[1])
lines <- lines[-1]

# ---- 2. Parse matrix ----
genomes <- character(n)
mat <- matrix(NA_real_, n, n)

for (i in seq_len(n)) {
  parts <- strsplit(lines[i], "\t")[[1]]
  genomes[i] <- parts[1]
  
  if (length(parts) > 1) {
    vals <- as.numeric(parts[-1])
    mat[i, 1:(i-1)] <- vals
  }
}

# ---- 3. Symmetrize ----

mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
diag(mat) <- 100

rownames(mat) <- basename(genomes)
colnames(mat) <- basename(genomes)

# ---- 4. CLEAN DATA (CRITICAL STEP) ----

# Convert 0 → NA (failed comparisons)
mat[mat == 0] <- NA

# Remove genomes with too many NA (broken ones)
keep <- rowSums(!is.na(mat)) > (0.5 * ncol(mat))
mat <- mat[keep, keep]

# ---- 5. Convert to distance ----
dist_mat <- 100 - mat

# Replace remaining NA with SLIGHTLY larger than max distance
max_val <- max(dist_mat, na.rm = TRUE)
dist_mat[is.na(dist_mat)] <- max_val + 0.01

# ---- 6. Add tiny noise (prevents identical values issue) ----
set.seed(1)
dist_mat <- dist_mat + matrix(rnorm(length(dist_mat), 0, 1e-6),
                              nrow(dist_mat))

# ---- 7. Plot ----
hmap <- pheatmap(
  mat,
  cluster_rows = TRUE,
  cellwidth = 5,
  cellheight = 5,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "orange", "red"))(100),
  main = "ANI Distance Heatmap",
  fontsize = 5
)

ggsave("hmap.svg", plot = hmap)

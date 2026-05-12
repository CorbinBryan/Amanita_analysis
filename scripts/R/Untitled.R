library(ggtree)
library(treeio)
library(ggplot2)
library(dplyr)
library(processx)
library(tidyverse)
library(ggtreeExtra)
library(ape)
setwd("/Users/corbinbryan/Desktop/Amanita_analysis/scripts/R")
treefile <- "../../astral4_anno2.tree"
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
  select(node, q1, q2, q3) #%>%
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

pies <- nodepie(freq_long, cols = c("q1", "q2", "q3"), 
                color = c("q1" = "#377eb8", "q2" = "#e41a1c", "q3" = "#4daf4a"))
p <- ggtree(tree_rooted) + 
  #geom_tiplab() + 
  theme_tree() + 
  geom_cladelabel(node=77, label='italic("A. muscaria")', parse = TRUE, color = "red") + 
  geom_cladelabel(node=101, label='italic("A. chrysoblema")',parse = TRUE, color = "skyblue") + 
  geom_cladelabel(node=99, label='italic("A. flavivolvata")',parse = TRUE, color = "gold") +
  xlim(0, max(tree_rooted@phylo[["edge.length"]]) * 1.1) 
j <- p %<+% name_key + geom_tiplab(aes(subset=(!label %in% exclude),
                                  label=paste0("italic('", species, "')"), 
                                  parse = TRUE), parse = TRUE) 


# Add pie charts as insets
final_plot <- inset(j, pies, width = 0.02, height = 0.02)


print(final_plot) 

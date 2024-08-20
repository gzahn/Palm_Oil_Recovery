# BUILD PHYLOQGENIES ####


# SETUP ####
set.seed(666)
mylib <- "/uufs/chpc.utah.edu/common/home/u6033249/R/library-4.4"

## packages ####
library(tidyverse, lib.loc = mylib)
library(phyloseq, lib.loc = mylib)
library(DECIPHER, lib.loc = mylib)
library(janitor, lib.loc = mylib)
library(phylogram, lib.loc = mylib) # on CRAN
# library(tidyverse)
# library(phyloseq)
# library(DECIPHER)
# library(phylogram)


## functions ####
source("./R/functions.R")

## load metadata ####
meta <- read_csv("./data/Metadata.csv") %>% 
  janitor::clean_names()

## load physeq objects ####
ps_its <- readRDS("./output/ITS_Physeq_cleaned.RDS")
ps_16s <- readRDS("./output/16S_Physeq_cleaned.RDS")

## extract sequences ####

seqs_16s <- rownames(tax_table(ps_16s))
names(seqs_16s) <- paste0("ASV_",1:length(seqs_16s)) # This propagates to the tip labels of the tree
seqs_its <- rownames(tax_table(ps_its))
names(seqs_its) <- paste0("ASV_",1:length(seqs_its)) # This propagates to the tip labels of the tree



# ALIGNMENT ####
alignment_16s <- DECIPHER::AlignSeqs(DNAStringSet(seqs_16s),processors=NULL)
alignment_its <- DECIPHER::AlignSeqs(DNAStringSet(seqs_its),processors=NULL)

# save progress 
saveRDS(alignment_16s,"./output/16S_dna_alignment.RDS")
alignment_16s <- readRDS("./Output/16S_dna_alignment_muscle.RDS")
saveRDS(alignment_its,"./output/ITS_dna_alignment.RDS")
alignment_its <- readRDS("./output/ITS_dna_alignment.RDS")

# DISTANCE MATRIX ####
dist_16s <- DistanceMatrix(alignment_16s,type = "dist",correction="Jukes-Cantor", verbose=FALSE,processors = NULL)
dist_its <- DistanceMatrix(alignment_its,type = "dist",correction="Jukes-Cantor", verbose=FALSE,processors = NULL)
# # distm <- DistanceMatrix(alignment,type = "matrix",correction="Jukes-Cantor", verbose=FALSE,processors = NULL)
# names(seqs) <- paste0("ASV_",1:length(seqs)) # This propagates to the tip labels of the tree

# BUILD INITIAL TREES ####
# Nieghbor-joining tree
treeNJ_16s <- TreeLine(myDistMatrix=dist_16s, method="NJ", cutoff=0.05, showPlot=FALSE, verbose=FALSE)
treeNJ_its <- TreeLine(myDistMatrix=dist_its, method="NJ", cutoff=0.05, showPlot=FALSE, verbose=FALSE)

# save progress
saveRDS(treeNJ_16s, "./output/16S_treeNJ.RDS")
treeNJ_16s <- readRDS("./output/16S_treeNJ.RDS")
saveRDS(treeNJ_its, "./output/ITS_treeNJ.RDS")
treeNJ_its <- readRDS("./output/ITS_treeNJ.RDS")

# convert to phylogram
tree_16s <- phylogram::as.phylo(treeNJ_16s)
tree_its <- phylogram::as.phylo(treeNJ_its)

# make sure tip labels match sample name order
asv_order_its <- tree_its$tip.label %>% str_remove("ASV_") %>% as.numeric
asv_order_16s <- tree_16s$tip.label %>% str_remove("ASV_") %>% as.numeric
tree_16s$tip.label <- taxa_names(ps_16s)[asv_order_16s]
tree_its$tip.label <- taxa_names(ps_its)[asv_order_its]

# ADD TREES TO PHYSEQ ####
# add trees to physeq objects
ps_16s_w_tree <- phyloseq(otu_table(ps_16s),
                          tax_table(ps_16s),
                          sample_data(ps_16s),
                          phy_tree(tree_16s))
ps_its_w_tree <- phyloseq(otu_table(ps_its),
                          tax_table(ps_its),
                          sample_data(ps_its),
                          phy_tree(tree_its))

# ps_16s_w_tree <- readRDS("./output/16S_Physeq_cleaned_w_tree.RDS")
# ps_its_w_tree <- readRDS("./output/ITS_Physeq_cleaned_w_tree.RDS")
# ps_16s_w_tree@sam_data <- ps_16s@sam_data
# ps_its_w_tree@sam_data <- ps_its@sam_data


# EXPORT ####
saveRDS(ps_16s_w_tree,"./output/16S_Physeq_cleaned_w_tree.RDS")
saveRDS(ps_its_w_tree,"./output/ITS_Physeq_cleaned_w_tree.RDS")


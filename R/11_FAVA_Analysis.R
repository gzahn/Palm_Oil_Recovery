# RUN FAVA ON TIME SERIES
# https://maikemorrison.github.io/FAVA/articles/microbiome_tutorial.html
# The FAVA R package implements the statistic FAVA, an FST-based Assessment of Variability across vectors of relative Abundances


# note: samples are nested within times!

# SETUP ####
set.seed(666)

## packages ####
library(tidyverse)
library(phyloseq)
library(FAVA)

## functions ####
source("./R/functions.R")

## load physeq objects ####
fung <- readRDS("./output/ITS_Physeq_cleaned_w_tree.RDS")
# import bacterial data glommed to species (not ASV level)
bact <- readRDS("./output/16S_Physeq_cleaned_w_tree.RDS") %>% tax_glom(taxrank = "Species",NArm = FALSE)

# FORMAT FOR FAVA ####

## generate metadata+relabund matrices ####
fung_dat <- FAVA::relab_phyloseq(fung)
bact_dat <- FAVA::relab_phyloseq(bact)

## build phylogenetic dissimilarity matrices ####
fung_tree <- phy_tree(fung)
bact_tree <- phy_tree(bact)

fung_dist <- ape::cophenetic.phylo(fung_tree)
bact_dist <- ape::cophenetic.phylo(bact_tree)

# Get the names of the species in your relative abundance matrix
nmetacols <- 1 + length(colnames(fung@sam_data))
fung_order <- colnames(fung_dat[,c(nmetacols:ncol(fung_dat))])
bact_order <- colnames(bact_dat[,c(nmetacols:ncol(bact_dat))])

# Confirm that the entries of the similarity matrix
# correspond to relative abundance matrix
if(!all(
  all(fung_order == colnames(fung_dist)),
  all(fung_order == rownames(fung_dist)),
  all(bact_order == colnames(bact_dist)),
  all(bact_order == rownames(bact_dist))
)){
  cat("Error: Check and correct taxa names and orders")
  }

## build similarity matrices ####
fung_simil_mat <- exp(-fung_dist)
bact_simil_mat <- exp(-bact_dist)

# are all diag elements still == 1 ?
diag(fung_simil_mat) %>% summary # eh?
diag(bact_simil_mat) %>% summary # eh?

# quick heatmap sanity check
heatmap(fung_simil_mat)

# VISUALIZE RELABUNDS ####
species_palette <- viridis::turbo(nrow(fung_simil_mat))[sample(1:nrow(fung_simil_mat))] %>%
  `names<-`(colnames(fung_dat)[nmetacols:ncol(fung_dat)])

# Make a ggplot2 stacked bar plot
plot_relabund(fung_dat,
              group = "treatment",
              arrange = "both",
              K = length(nmetacols:ncol(fung_dat))) +
  # Specify a custom color scheme
  ggplot2::scale_color_manual(values = species_palette) +
  ggplot2::scale_fill_manual(values = species_palette)

# COMPUTE FAVA STAT ####
# unweighted
fung_fava_uw <- 
  fava(relab_matrix = fung_dat,
     group = "treatment",
     K = length(nmetacols:ncol(fung_dat)))
bact_fava_uw <- 
  fava(relab_matrix = bact_dat,
       group = "treatment",
       K = length(nmetacols:ncol(bact_dat)))
# weighted
fung_fava_w <- 
  fava(relab_matrix = fung_dat,
     group = "treatment",
     K = length(nmetacols:ncol(fung_dat)),
     S = fung_simil_mat)
bact_fava_w <- 
  fava(relab_matrix = bact_dat,
       group = "treatment",
       K = length(nmetacols:ncol(bact_dat)),
       S = bact_simil_mat)

# BOOTSTRAPPING ####
# weighted by phlogenetic similarity
fung_fava_bs <- bootstrap_fava(relab_matrix = fung_dat,
                               n_replicates = 1000,
                               seed = 666,
                               group = "treatment",
                               K = length(nmetacols:ncol(fung_dat)),
                               S = fung_simil_mat)
bact_fava_bs <- bootstrap_fava(relab_matrix = bact_dat,
                               n_replicates = 1000,
                               seed = 666,
                               group = "treatment",
                               K = length(nmetacols:ncol(bact_dat)),
                               S = bact_simil_mat)

# save FAVA bootstrap output, since it takes forever to run
saveRDS(fung_fava_bs,"./output/bootstrap_fava_fungi.RDS")
saveRDS(bact_fava_bs,"./output/bootstrap_fava_bacteria.RDS")

# examine p-values
fung_fava_bs$P_values %>% 
  full_join(fung_fava_bs$observed_difference)
bact_fava_bs$P_values



# EXPORT RESULTS AND FIGS ####

# FAVA stats tables
fung_fava_bs$observed_difference %>% 
  mutate(Comparison = Comparison %>% str_replace("\n"," ")) %>% 
  full_join(fung_fava_bs$P_values) %>% 
  write_csv("./output/fava_stats_fungi.csv")
bact_fava_bs$observed_difference %>% 
  mutate(Comparison = Comparison %>% str_replace("\n"," ")) %>% 
  full_join(bact_fava_bs$P_values) %>% 
  write_csv("./output/fava_stats_bacteria.csv")

# plots
ggsave(plot = fung_fava_bs$bootstrap_distribution_plot,
       filename = "./output/figs/fava_bootstrap_plot_fungi.png",
       dpi=400,
       width = 16,
       height = 8)
ggsave(plot = bact_fava_bs$bootstrap_distribution_plot,
       filename = "./output/figs/fava_bootstrap_plot_bacteria.png",
       dpi=400,
       width = 16,
       height = 8)

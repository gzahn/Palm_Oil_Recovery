# NETWORK ANALYSIS

# SETUP ####

## packages ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(broom); packageVersion("broom")
library(patchwork); packageVersion("patchwork")
library(igraph); packageVersion("igraph")
library(SpiecEasi); packageVersion("SpiecEasi")
library(hubfindr); packageVersion("hubfindr") # devtools::install_github("gzahn/hubfindr")
library(corncob); packageVersion("corncob")

## functions
source("./R/functions.R")

## options ####
options(warn=0)
set.seed(666)

theme_set(theme_bw() +
            theme(axis.title = element_text(face='bold',size=18),
                  axis.text = element_text(face='bold',size=14),
                  legend.title = element_text(face='bold',size=18),
                  legend.text = element_text(face='bold',size=14),
                  strip.text = element_text(face='bold',size=18),
                  plot.title = element_text(face='bold',size=18)))

## load physeq objects ####
# agglomerated taxa to species
fung <- readRDS("./output/ITS_Physeq_cleaned_w_tree.RDS") %>% tax_glom("Species",NArm = FALSE)
bact <- readRDS("./output/16S_Physeq_cleaned_w_tree.RDS") %>% tax_glom("Species",NArm = FALSE)

# full physeq with fungi and bacteria
# make copies
fung2 <- fung
bact2 <- bact
# remove phy trees
fung2@phy_tree <- NULL
bact2@phy_tree <- NULL
# make sample names match
sample_names(fung2) <- sample_names(fung2) %>% str_remove("ITS_")
sample_names(bact2) <- sample_names(bact2) %>% str_remove("16S_")
# merge to single physeq
full <- merge_phyloseq(fung2,bact2)


# BUILD NETWORKS ####
# (all networks built on species-level consolidations of ASVs)

hubfindr::find_hubs

# SpiecEasi parameters
se.params <- list(rep.num=20, ncores=(parallel::detectCores()-1))

## Fungi ####
se.mb.fung <- SpiecEasi::spiec.easi(data = fung,
                                    method='mb',
                                    sel.criterion = "bstars",
                                    pulsar.params=se.params)
saveRDS(se.mb.fung,"./output/Fungal_SpiecEasi_out.RDS")
se.mb.fung <- readRDS("./output/Fungal_SpiecEasi_out.RDS")

### per-treatment × domain networks ####
# oil palm
fung_oil_palm <- subset_samples(fung, treatment == "oil palm")
fung_oil_palm <- prune_taxa(taxa_sums(fung_oil_palm) > 0, fung_oil_palm)
se.mb.fung_oil_palm <- SpiecEasi::spiec.easi(data = fung_oil_palm, method='mb',
                                             sel.criterion = "bstars", pulsar.params=se.params)
saveRDS(se.mb.fung_oil_palm, "./output/Fungal_SpiecEasi_oil_palm.RDS")

# 0-1 years
fung_0_1yrs <- subset_samples(fung, treatment == "0-1 years")
fung_0_1yrs <- prune_taxa(taxa_sums(fung_0_1yrs) > 0, fung_0_1yrs)
se.mb.fung_0_1yrs <- SpiecEasi::spiec.easi(data = fung_0_1yrs, method='mb',
                                           sel.criterion = "bstars", pulsar.params=se.params)
saveRDS(se.mb.fung_0_1yrs, "./output/Fungal_SpiecEasi_0_1yrs.RDS")

# 1-2 years
fung_1_2yrs <- subset_samples(fung, treatment == "1-2 years")
fung_1_2yrs <- prune_taxa(taxa_sums(fung_1_2yrs) > 0, fung_1_2yrs)
se.mb.fung_1_2yrs <- SpiecEasi::spiec.easi(data = fung_1_2yrs, method='mb',
                                           sel.criterion = "bstars", pulsar.params=se.params)
saveRDS(se.mb.fung_1_2yrs, "./output/Fungal_SpiecEasi_1_2yrs.RDS")

# 2-3 years
fung_2_3yrs <- subset_samples(fung, treatment == "2-3 years")
fung_2_3yrs <- prune_taxa(taxa_sums(fung_2_3yrs) > 0, fung_2_3yrs)
se.mb.fung_2_3yrs <- SpiecEasi::spiec.easi(data = fung_2_3yrs, method='mb',
                                           sel.criterion = "bstars", pulsar.params=se.params)
saveRDS(se.mb.fung_2_3yrs, "./output/Fungal_SpiecEasi_2_3yrs.RDS")

# 3-4 years
fung_3_4yrs <- subset_samples(fung, treatment == "3-4 years")
fung_3_4yrs <- prune_taxa(taxa_sums(fung_3_4yrs) > 0, fung_3_4yrs)
se.mb.fung_3_4yrs <- SpiecEasi::spiec.easi(data = fung_3_4yrs, method='mb',
                                           sel.criterion = "bstars", pulsar.params=se.params)
saveRDS(se.mb.fung_3_4yrs, "./output/Fungal_SpiecEasi_3_4yrs.RDS")


## Bacteria ####
se.mb.bact <- SpiecEasi::spiec.easi(data = bact,
                                    method='mb',
                                    sel.criterion = "bstars",
                                    pulsar.params=se.params)
saveRDS(se.mb.bact,"./output/Bacterial_SpiecEasi_out.RDS")
se.mb.bact <- readRDS("./output/Bacterial_SpiecEasi_out.RDS")

### per-treatment × domain networks ####
# oil palm
bact_oil_palm <- subset_samples(bact, treatment == "oil palm")
bact_oil_palm <- prune_taxa(taxa_sums(bact_oil_palm) > 0, bact_oil_palm)
se.mb.bact_oil_palm <- SpiecEasi::spiec.easi(data = bact_oil_palm, method = 'mb',
                                             sel.criterion = "bstars", pulsar.params = se.params)
saveRDS(se.mb.bact_oil_palm, "./output/Bacterial_SpiecEasi_oil_palm.RDS")

# 0-1 years
bact_0_1yrs <- subset_samples(bact, treatment == "0-1 years")
bact_0_1yrs <- prune_taxa(taxa_sums(bact_0_1yrs) > 0, bact_0_1yrs)
se.mb.bact_0_1yrs <- SpiecEasi::spiec.easi(data = bact_0_1yrs, method = 'mb',
                                           sel.criterion = "bstars", pulsar.params = se.params)
saveRDS(se.mb.bact_0_1yrs, "./output/Bacterial_SpiecEasi_0_1yrs.RDS")

# 1-2 years
bact_1_2yrs <- subset_samples(bact, treatment == "1-2 years")
bact_1_2yrs <- prune_taxa(taxa_sums(bact_1_2yrs) > 0, bact_1_2yrs)
se.mb.bact_1_2yrs <- SpiecEasi::spiec.easi(data = bact_1_2yrs, method = 'mb',
                                           sel.criterion = "bstars", pulsar.params = se.params)
saveRDS(se.mb.bact_1_2yrs, "./output/Bacterial_SpiecEasi_1_2yrs.RDS")

# 2-3 years
bact_2_3yrs <- subset_samples(bact, treatment == "2-3 years")
bact_2_3yrs <- prune_taxa(taxa_sums(bact_2_3yrs) > 0, bact_2_3yrs)
se.mb.bact_2_3yrs <- SpiecEasi::spiec.easi(data = bact_2_3yrs, method = 'mb',
                                           sel.criterion = "bstars", pulsar.params = se.params)
saveRDS(se.mb.bact_2_3yrs, "./output/Bacterial_SpiecEasi_2_3yrs.RDS")

# 3-4 years
bact_3_4yrs <- subset_samples(bact, treatment == "3-4 years")
bact_3_4yrs <- prune_taxa(taxa_sums(bact_3_4yrs) > 0, bact_3_4yrs)
se.mb.bact_3_4yrs <- SpiecEasi::spiec.easi(data = bact_3_4yrs, method = 'mb',
                                           sel.criterion = "bstars", pulsar.params = se.params)
saveRDS(se.mb.bact_3_4yrs, "./output/Bacterial_SpiecEasi_3_4yrs.RDS")


## Both ####
se.mb.full <- SpiecEasi::spiec.easi(data = full,
                                    method='mb',
                                    sel.criterion = "bstars",
                                    pulsar.params=se.params)
saveRDS(se.mb.full,"./output/Full_SpiecEasi_out.RDS")
se.mb.full <- readRDS("./output/Full_SpiecEasi_out.RDS")

### per-treatment × domain networks ####
# oil palm
full_oil_palm <- subset_samples(full, treatment == "oil palm")
full_oil_palm <- prune_taxa(taxa_sums(full_oil_palm) > 0, full_oil_palm)
se.mb.full_oil_palm <- SpiecEasi::spiec.easi(data = full_oil_palm, method = 'mb',
                                             sel.criterion = "bstars", pulsar.params = se.params)
saveRDS(se.mb.full_oil_palm, "./output/Full_SpiecEasi_oil_palm.RDS")

# 0-1 years
full_0_1yrs <- subset_samples(full, treatment == "0-1 years")
full_0_1yrs <- prune_taxa(taxa_sums(full_0_1yrs) > 0, full_0_1yrs)
se.mb.full_0_1yrs <- SpiecEasi::spiec.easi(data = full_0_1yrs, method = 'mb',
                                           sel.criterion = "bstars", pulsar.params = se.params)
saveRDS(se.mb.full_0_1yrs, "./output/Full_SpiecEasi_0_1yrs.RDS")

# 1-2 years
full_1_2yrs <- subset_samples(full, treatment == "1-2 years")
full_1_2yrs <- prune_taxa(taxa_sums(full_1_2yrs) > 0, full_1_2yrs)
se.mb.full_1_2yrs <- SpiecEasi::spiec.easi(data = full_1_2yrs, method = 'mb',
                                           sel.criterion = "bstars", pulsar.params = se.params)
saveRDS(se.mb.full_1_2yrs, "./output/Full_SpiecEasi_1_2yrs.RDS")

# 2-3 years
full_2_3yrs <- subset_samples(full, treatment == "2-3 years")
full_2_3yrs <- prune_taxa(taxa_sums(full_2_3yrs) > 0, full_2_3yrs)
se.mb.full_2_3yrs <- SpiecEasi::spiec.easi(data = full_2_3yrs, method = 'mb',
                                           sel.criterion = "bstars", pulsar.params = se.params)
saveRDS(se.mb.full_2_3yrs, "./output/Full_SpiecEasi_2_3yrs.RDS")

# 3-4 years
full_3_4yrs <- subset_samples(full, treatment == "3-4 years")
full_3_4yrs <- prune_taxa(taxa_sums(full_3_4yrs) > 0, full_3_4yrs)
se.mb.full_3_4yrs <- SpiecEasi::spiec.easi(data = full_3_4yrs, method = 'mb',
                                           sel.criterion = "bstars", pulsar.params = se.params)
saveRDS(se.mb.full_3_4yrs, "./output/Full_SpiecEasi_3_4yrs.RDS")

## build igraph ####
fung_igraph <- adj2igraph(getRefit(se.mb.fung),vertex.attr = list(name=taxa_names(fung)))
bact_igraph <- adj2igraph(getRefit(se.mb.bact),vertex.attr = list(name=taxa_names(bact)))
full_igraph <- adj2igraph(getRefit(se.mb.full),vertex.attr = list(name=taxa_names(full)))

### for treatment-level graphs ####
fung_igraph_oilpalm <- adj2igraph(getRefit(se.mb.fung_oil_palm),vertex.attr = list(name=taxa_names(fung_oil_palm)))
fung_igraph_01 <- adj2igraph(getRefit(se.mb.fung_0_1yrs),vertex.attr = list(name=taxa_names(fung_0_1yrs)))
fung_igraph_12 <- adj2igraph(getRefit(se.mb.fung_1_2yrs),vertex.attr = list(name=taxa_names(fung_1_2yrs)))
fung_igraph_23 <- adj2igraph(getRefit(se.mb.fung_2_3yrs),vertex.attr = list(name=taxa_names(fung_2_3yrs)))
fung_igraph_34 <- adj2igraph(getRefit(se.mb.fung_3_4yrs),vertex.attr = list(name=taxa_names(fung_3_4yrs)))

bact_igraph_oilpalm <- adj2igraph(getRefit(se.mb.bact_oil_palm),vertex.attr = list(name=taxa_names(bact_oil_palm)))
bact_igraph_01 <- adj2igraph(getRefit(se.mb.bact_0_1yrs),vertex.attr = list(name=taxa_names(bact_0_1yrs)))
bact_igraph_12 <- adj2igraph(getRefit(se.mb.bact_1_2yrs),vertex.attr = list(name=taxa_names(bact_1_2yrs)))
bact_igraph_23 <- adj2igraph(getRefit(se.mb.bact_2_3yrs),vertex.attr = list(name=taxa_names(bact_2_3yrs)))
bact_igraph_34 <- adj2igraph(getRefit(se.mb.bact_3_4yrs),vertex.attr = list(name=taxa_names(bact_3_4yrs)))

full_igraph_oilpalm <- adj2igraph(getRefit(se.mb.full_oil_palm),vertex.attr = list(name=taxa_names(full_oil_palm)))
full_igraph_01 <- adj2igraph(getRefit(se.mb.full_0_1yrs),vertex.attr = list(name=taxa_names(full_0_1yrs)))
full_igraph_12 <- adj2igraph(getRefit(se.mb.full_1_2yrs),vertex.attr = list(name=taxa_names(full_1_2yrs)))
full_igraph_23 <- adj2igraph(getRefit(se.mb.full_2_3yrs),vertex.attr = list(name=taxa_names(full_2_3yrs)))
full_igraph_34 <- adj2igraph(getRefit(se.mb.full_3_4yrs),vertex.attr = list(name=taxa_names(full_3_4yrs)))

# Save global igraphs
saveRDS(fung_igraph, "./output/fungal_igraph_full.RDS")
saveRDS(bact_igraph, "./output/bacterial_igraph_full.RDS")
saveRDS(full_igraph, "./output/combined_igraph_full.RDS")

# Save treatment-level fungal igraphs
saveRDS(fung_igraph_oilpalm, "./output/fungal_igraph_oil_palm.RDS")
saveRDS(fung_igraph_01,       "./output/fungal_igraph_0_1yrs.RDS")
saveRDS(fung_igraph_12,       "./output/fungal_igraph_1_2yrs.RDS")
saveRDS(fung_igraph_23,       "./output/fungal_igraph_2_3yrs.RDS")
saveRDS(fung_igraph_34,       "./output/fungal_igraph_3_4yrs.RDS")

# Save treatment-level bacterial igraphs
saveRDS(bact_igraph_oilpalm, "./output/bacterial_igraph_oil_palm.RDS")
saveRDS(bact_igraph_01,      "./output/bacterial_igraph_0_1yrs.RDS")
saveRDS(bact_igraph_12,      "./output/bacterial_igraph_1_2yrs.RDS")
saveRDS(bact_igraph_23,      "./output/bacterial_igraph_2_3yrs.RDS")
saveRDS(bact_igraph_34,      "./output/bacterial_igraph_3_4yrs.RDS")

# Save treatment-level combined igraphs
saveRDS(full_igraph_oilpalm, "./output/combined_igraph_oil_palm.RDS")
saveRDS(full_igraph_01,      "./output/combined_igraph_0_1yrs.RDS")
saveRDS(full_igraph_12,      "./output/combined_igraph_1_2yrs.RDS")
saveRDS(full_igraph_23,      "./output/combined_igraph_2_3yrs.RDS")
saveRDS(full_igraph_34,      "./output/combined_igraph_3_4yrs.RDS")

# save subset physeq objects
# Save fungal phyloseq objects
saveRDS(fung_oil_palm, "./output/fungal_phyloseq_oil_palm.RDS")
saveRDS(fung_0_1yrs,   "./output/fungal_phyloseq_0_1yrs.RDS")
saveRDS(fung_1_2yrs,   "./output/fungal_phyloseq_1_2yrs.RDS")
saveRDS(fung_2_3yrs,   "./output/fungal_phyloseq_2_3yrs.RDS")
saveRDS(fung_3_4yrs,   "./output/fungal_phyloseq_3_4yrs.RDS")

# Save bacterial phyloseq objects
saveRDS(bact_oil_palm, "./output/bacterial_phyloseq_oil_palm.RDS")
saveRDS(bact_0_1yrs,   "./output/bacterial_phyloseq_0_1yrs.RDS")
saveRDS(bact_1_2yrs,   "./output/bacterial_phyloseq_1_2yrs.RDS")
saveRDS(bact_2_3yrs,   "./output/bacterial_phyloseq_2_3yrs.RDS")
saveRDS(bact_3_4yrs,   "./output/bacterial_phyloseq_3_4yrs.RDS")

# Save combined (fungi + bacteria) phyloseq objects
saveRDS(full_oil_palm, "./output/combined_phyloseq_oil_palm.RDS")
saveRDS(full_0_1yrs,   "./output/combined_phyloseq_0_1yrs.RDS")
saveRDS(full_1_2yrs,   "./output/combined_phyloseq_1_2yrs.RDS")
saveRDS(full_2_3yrs,   "./output/combined_phyloseq_2_3yrs.RDS")
saveRDS(full_3_4yrs,   "./output/combined_phyloseq_3_4yrs.RDS")

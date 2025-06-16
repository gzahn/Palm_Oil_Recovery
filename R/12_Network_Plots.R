# NETWORK PLOTS 

# SETUP ####

## packages ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(broom); packageVersion("broom")
library(patchwork); packageVersion("patchwork")
library(igraph); packageVersion("igraph")

source("./R/functions.R")

# LOAD DATA ####

## Load physeq + traits ####

# Load cleaned phyloseq object
bact <- readRDS("./output/16S_Physeq_cleaned_w_tree.RDS")
fung <- readRDS("./output/ITS_Physeq_cleaned_w_tree.RDS")

# load trait data frames
fung_trait <- readRDS("./output/fungal_trait_df.RDS")
bact_trait <- readRDS("./output/bacterial_enzyme_diversity_df.RDS")

# add known enzyme diversity to kingdom slot of bact
bact_genus <- 
bact %>% 
  tax_table() %>% 
  as.data.frame() %>% 
  pluck('Genus')

bact@tax_table[,1] <- 
data.frame(genus=bact_genus) %>% 
  left_join(bact_trait) %>% 
  pluck("num_known_enzymes") %>% 
  as.character()

# add fungal trophic guild to fung kingdom slot
fung_asv <- 
  fung %>% 
  taxa_names()
fung@tax_table[,1] <- 
  data.frame(ASV=fung_asv) %>% 
  left_join(fung_trait) %>% 
  pluck("trophicMode") %>% 
  as.character()



## Load network objects ####

# Read treatment-level fungal igraphs
fung_igraph_oilpalm <- readRDS("./output/fungal_igraph_oil_palm.RDS")
fung_igraph_01 <- readRDS("./output/fungal_igraph_0_1yrs.RDS")
fung_igraph_12 <- readRDS("./output/fungal_igraph_1_2yrs.RDS")
fung_igraph_23 <- readRDS("./output/fungal_igraph_2_3yrs.RDS")
fung_igraph_34 <- readRDS("./output/fungal_igraph_3_4yrs.RDS")

# Read treatment-level bacterial igraphs
bact_igraph_oilpalm <- readRDS("./output/bacterial_igraph_oil_palm.RDS")
bact_igraph_01 <- readRDS("./output/bacterial_igraph_0_1yrs.RDS")
bact_igraph_12 <- readRDS("./output/bacterial_igraph_1_2yrs.RDS")
bact_igraph_23 <- readRDS("./output/bacterial_igraph_2_3yrs.RDS")
bact_igraph_34 <- readRDS("./output/bacterial_igraph_3_4yrs.RDS")

# Read treatment-level combined igraphs
full_igraph_oilpalm <- readRDS("./output/combined_igraph_oil_palm.RDS")
full_igraph_01 <- readRDS("./output/combined_igraph_0_1yrs.RDS")
full_igraph_12 <- readRDS("./output/combined_igraph_1_2yrs.RDS")
full_igraph_23 <- readRDS("./output/combined_igraph_2_3yrs.RDS")
full_igraph_34 <- readRDS("./output/combined_igraph_3_4yrs.RDS")



# ADD TRAITS TO igraph ####


for(i in ls(pattern = "full_igraph_")){
  x <- get(i)
  print(i)
  # get vertice names
  graph.asvs <- vertex.attributes(x)[['name']]
  # find which are bacterial
  graph.bact.vertices <- graph.asvs %in% taxa_names(bact)
  # color based on bact vs fungi
  V(x)$color <- ifelse(graph.bact.vertices,pal$pal.earthtones[1],pal$pal.earthtones[2])
  # plot
  plot_hubs(x,layout = 'auto')
}

for(i in ls(pattern = "fung_igraph_")){
  x <- get(i)
  print(i)
  # get vertice names
  graph.asvs <- vertex.attributes(x)[['name']]
  # find which are bacterial
  graph.bact.vertices <- graph.asvs %in% taxa_names(bact)
  # color based on bact vs fungi
  V(x)$color <- ifelse(graph.bact.vertices,pal$pal.earthtones[1],pal$pal.earthtones[2])
  # plot
  plot_hubs(x,layout = 'auto')
}

for(i in ls(pattern = "bact_igraph_")){
  x <- get(i)
  print(i)
  # get vertice names
  graph.asvs <- vertex.attributes(x)[['name']]
  # find which are bacterial
  graph.bact.vertices <- graph.asvs %in% taxa_names(bact)
  # color based on bact vs fungi
  V(x)$color <- ifelse(graph.bact.vertices,pal$pal.earthtones[1],pal$pal.earthtones[2])
  # plot
  plot_hubs(x,layout = 'auto')
}







# PLOT HUBS ####

# colored by trait

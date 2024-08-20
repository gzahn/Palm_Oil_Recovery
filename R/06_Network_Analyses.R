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

# SpiecEasi parameters
se.params <- list(rep.num=20, ncores=(parallel::detectCores()-1))

## Fungi ####
se.mb.fung <- SpiecEasi::spiec.easi(data = fung,
                                    method='mb',
                                    sel.criterion = "bstars",
                                    pulsar.params=se.params)
saveRDS(se.mb.fung,"./output/Fungal_SpiecEasi_out.RDS")
se.mb.fung <- readRDS("./output/Fungal_SpiecEasi_out.RDS")

## Bacteria ####
se.mb.bact <- SpiecEasi::spiec.easi(data = bact,
                                    method='mb',
                                    sel.criterion = "bstars",
                                    pulsar.params=se.params)
saveRDS(se.mb.fung,"./output/Bacterial_SpiecEasi_out.RDS")
se.mb.fung <- readRDS("./output/Bacterial_SpiecEasi_out.RDS")

## Both ####
se.mb.full <- SpiecEasi::spiec.easi(data = full,
                                    method='mb',
                                    sel.criterion = "bstars",
                                    pulsar.params=se.params)
saveRDS(se.mb.full,"./output/Full_SpiecEasi_out.RDS")
se.mb.full <- readRDS("./output/Full_SpiecEasi_out.RDS")


## build igraph ####
fung_igraph <- adj2igraph(getRefit(se.mb.fung),vertex.attr = list(name=taxa_names(fung)))
bact_igraph <- adj2igraph(getRefit(se.mb.bact),vertex.attr = list(name=taxa_names(bact)))
full_igraph <- adj2igraph(getRefit(se.mb.full),vertex.attr = list(name=taxa_names(full)))


# FIND HUB TAXA ####
fung_hub_taxa <- hubfindr::find_hubs(graph = fung_igraph,physeq = fung)
bact_hub_taxa <- hubfindr::find_hubs(graph = bact_igraph,physeq = bact)
full_hub_taxa <- hubfindr::find_hubs(graph = full_igraph,physeq = full)

## add taxonomy ####
fung_hub_taxa$hub_taxon <- 
  corncob::otu_to_taxonomy(fung_hub_taxa$vertex_id,fung)
bact_hub_taxa$hub_taxon <- 
  corncob::otu_to_taxonomy(bact_hub_taxa$vertex_id,bact)
full_hub_taxa$hub_taxon <- 
  corncob::otu_to_taxonomy(full_hub_taxa$vertex_id,full)
## export ####
full_join(
  fung_hub_taxa %>% 
    as.data.frame() %>%
    mutate(taxa = "Fungi"),
  bact_hub_taxa %>% 
    as.data.frame() %>%
    mutate(taxa = "Bacteria")
) %>% 
  full_join(full_hub_taxa %>% 
              as.data.frame() %>%
              mutate(taxa = "Both")
            ) %>% 
  write_csv("./output/hub_taxa_list.csv")


# NETWORK PROPERTIES ####

## subgraphs for each treatment ####

### fungi ####
fungi_netstats <- list()
for(i in (fung@sam_data$sample_id)){
  # build subset on the fly
  ss <- fung %>% subset_samples(sample_id == i)
  ss <- ss %>% subset_taxa(taxa_sums(ss) > 0)
  # find network stats of that subset
  stat_list <- find_ig_subset_attr(ps.subset = ss,ig.full = fung_igraph)
  fungi_netstats[[i]] <- stat_list
}

### bacteria ####
bacteria_netstats <- list()
for(i in (bact@sam_data$sample_id)){
  # build subset on the fly
  ss <- bact %>% subset_samples(sample_id == i)
  ss <- ss %>% subset_taxa(taxa_sums(ss) > 0)
  # find network stats of that subset
  stat_list <- find_ig_subset_attr(ps.subset = ss,ig.full = bact_igraph)
  bacteria_netstats[[i]] <- stat_list
}

### both ####
full_netstats <- list()
for(i in (full@sam_data$sample_id)){
  # build subset on the fly
  ss <- full %>% subset_samples(sample_id == i)
  ss <- ss %>% subset_taxa(taxa_sums(ss) > 0)
  # find network stats of that subset
  stat_list <- find_ig_subset_attr(ps.subset = ss,ig.full = full_igraph)
  full_netstats[[i]] <- stat_list
}

## export network statistics ####
saveRDS(fungi_netstats,"./output/fungi_network_stats_list.RDS")
saveRDS(bacteria_netstats,"./output/bacterial_network_stats_list.RDS")
saveRDS(full_netstats,"./output/full_network_stats_list.RDS")

## combine all network stats ####
# add all to single object
complete_netstats <- 
  list(fungi=fungi_netstats,
       bacteria=bacteria_netstats,
       full=full_netstats)

# nested list...access as follows:
complete_netstats$fungi$P1.1$n_vertices

## Extract attributes ####

### Degree Distributions ####
# isolate degree_distributions From all subsets
fung_degdist <- 
  data.frame(taxa = "Fungi",
             deg_dist = I(
               complete_netstats$fungi %>% 
                 map('deg_dist'))
  )
bact_degdist <- 
  data.frame(taxa = "Bacteria",
             deg_dist = I(
               complete_netstats$bacteria %>% 
                 map('deg_dist'))
  )
full_degdist <- 
  data.frame(taxa = "Both",
             deg_dist = I(
               complete_netstats$full %>% 
                 map('deg_dist'))
  )

deg_dist_df <- 
  fung_degdist %>% 
  full_join(bact_degdist) %>% 
  full_join(full_degdist)
# add grouping variable for each taxon
deg_dist_df$treatment <- rep((fung@sam_data$sample_id),3)

# make deg_dist lengths match
m <-
  deg_dist_df$deg_dist %>% 
  map_dbl(length) %>% 
  max
deg_dist_df$deg_dist <- lapply(deg_dist_df$deg_dist, "length<-",m)  

# make deg_dist lengths match
m <-
  deg_dist_df$deg_dist %>% 
  map_dbl(length) %>% 
  max
deg_dist_df$deg_dist <- lapply(deg_dist_df$deg_dist, "length<-",m)  
# ensure good names for list elements
names(deg_dist_df$deg_dist) <- paste0(c(rep('fung_',80),rep('bact_',80),rep('full_',80)),
                                      rep((fung@sam_data$sample_id),3))

deg_dist_df <- 
  deg_dist_df$deg_dist %>% 
  as.data.frame() %>% 
  pivot_longer(everything(),names_to = "sample_id") %>% 
  dplyr::filter(!is.na(value)) %>%
  mutate(taxa = sample_id %>% str_split("_") %>% map_chr(1),
         sample_id = sample_id %>% str_split("_") %>% map_chr(2))
# make factor for taxa
deg_dist_df$taxa <- 
  deg_dist_df$taxa %>% 
  str_replace("fung","Fungi") %>% 
  str_replace("bact","Bacteria") %>% 
  str_replace("full","Both") %>% 
  factor(levels=c("Bacteria","Fungi","Both"))
# quick plot
deg_dist_df %>% 
  ggplot(aes(x=value,fill=taxa)) +
  geom_density(alpha=.25) +
  facet_wrap(~taxa,scales='free') +
  scale_fill_manual(values = pal$pal.earthtones)

### max degree ####
fung_maxdeg <- 
  data.frame(taxa = "Fungi",
             max_degree = I(
               complete_netstats$fungi %>% 
                 map('max_degree'))
  ) %>% mutate(sample_id = row.names(.))
bact_maxdeg <- 
  data.frame(taxa = "Bacteria",
             max_degree = I(
               complete_netstats$bacteria %>% 
                 map('max_degree'))
  ) %>% mutate(sample_id = row.names(.))
full_maxdeg <- 
  data.frame(taxa = "Both",
             max_degree = I(
               complete_netstats$full %>% 
                 map('max_degree'))
  ) %>% mutate(sample_id = row.names(.))

# create dataframe
max_deg_df <- 
fung_maxdeg %>% 
  full_join(bact_maxdeg) %>% 
  full_join(full_maxdeg)

### N vertices ####
fung_nvertices <- 
  data.frame(taxa = "Fungi",
             n_vertices = I(
               complete_netstats$fungi %>% 
                 map('n_vertices'))
  ) %>% mutate(sample_id = row.names(.))
bact_nvertices <- 
  data.frame(taxa = "Bacteria",
             n_vertices = I(
               complete_netstats$bacteria %>% 
                 map('n_vertices'))
  ) %>% mutate(sample_id = row.names(.))
full_nvertices <- 
  data.frame(taxa = "Both",
             n_vertices = I(
               complete_netstats$full %>% 
                 map('n_vertices'))
  ) %>% mutate(sample_id = row.names(.))
n_vertices_df <- 
  fung_nvertices %>% 
  full_join(bact_nvertices) %>% 
  full_join(full_nvertices)

### mean dist ####
fung_mean_dist <- 
  data.frame(taxa = "Fungi",
             mean_dist = I(
               complete_netstats$fungi %>% 
                 map('mean_dist'))
  ) %>% mutate(sample_id = row.names(.))
bact_mean_dist <- 
  data.frame(taxa = "Bacteria",
             mean_dist = I(
               complete_netstats$bacteria %>% 
                 map('mean_dist'))
  ) %>% mutate(sample_id = row.names(.))
full_mean_dist <- 
  data.frame(taxa = "Both",
             mean_dist = I(
               complete_netstats$full %>% 
                 map('mean_dist'))
  ) %>% mutate(sample_id = row.names(.))
mean_dist_df <- 
  fung_mean_dist %>% 
  full_join(bact_mean_dist) %>% 
  full_join(full_mean_dist)

### clique number ####
fung_clique_num <- 
  data.frame(taxa = "Fungi",
             clique_num = I(
               complete_netstats$fungi %>% 
                 map('clique_num'))
  ) %>% mutate(sample_id = row.names(.))
bact_clique_num <- 
  data.frame(taxa = "Bacteria",
             clique_num = I(
               complete_netstats$bacteria %>% 
                 map('clique_num'))
  ) %>% mutate(sample_id = row.names(.))
full_clique_num <- 
  data.frame(taxa = "Both",
             clique_num = I(
               complete_netstats$full %>% 
                 map('clique_num'))
  ) %>% mutate(sample_id = row.names(.))
clique_num_df <- 
  fung_clique_num %>% 
  full_join(bact_clique_num) %>% 
  full_join(full_clique_num)

### mean betweenness ####
fung_mean_betweenness <- 
  data.frame(taxa = "Fungi",
             mean_betweenness = I(
               complete_netstats$fungi %>% 
                 map('mean_betweenness'))
  ) %>% mutate(sample_id = row.names(.))
bact_mean_betweenness <- 
  data.frame(taxa = "Bacteria",
             mean_betweenness = I(
               complete_netstats$bacteria %>% 
                 map('mean_betweenness'))
  ) %>% mutate(sample_id = row.names(.))
full_mean_betweenness <- 
  data.frame(taxa = "Both",
             mean_betweenness = I(
               complete_netstats$full %>% 
                 map('mean_betweenness'))
  ) %>% mutate(sample_id = row.names(.))
mean_betweenness_df <- 
  fung_mean_betweenness %>% 
  full_join(bact_mean_betweenness) %>% 
  full_join(full_mean_betweenness)

### mean closeness ####
fung_mean_closeness <- 
  data.frame(taxa = "Fungi",
             mean_closeness = I(
               complete_netstats$fungi %>% 
                 map('mean_closeness'))
  ) %>% mutate(sample_id = row.names(.))
bact_mean_closeness <- 
  data.frame(taxa = "Bacteria",
             mean_closeness = I(
               complete_netstats$bacteria %>% 
                 map('mean_closeness'))
  ) %>% mutate(sample_id = row.names(.))
full_mean_closeness <- 
  data.frame(taxa = "Both",
             mean_closeness = I(
               complete_netstats$full %>% 
                 map('mean_closeness'))
  ) %>% mutate(sample_id = row.names(.))
mean_closeness_df <- 
  fung_mean_closeness %>% 
  full_join(bact_mean_closeness) %>% 
  full_join(full_mean_closeness)

### mean coreness ####
fung_mean_coreness <- 
  data.frame(taxa = "Fungi",
             mean_coreness = I(
               complete_netstats$fungi %>% 
                 map('mean_coreness'))
  ) %>% mutate(sample_id = row.names(.))
bact_mean_coreness <- 
  data.frame(taxa = "Bacteria",
             mean_coreness = I(
               complete_netstats$bacteria %>% 
                 map('mean_coreness'))
  ) %>% mutate(sample_id = row.names(.))
full_mean_coreness <- 
  data.frame(taxa = "Both",
             mean_coreness = I(
               complete_netstats$full %>% 
                 map('mean_coreness'))
  ) %>% mutate(sample_id = row.names(.))
mean_coreness_df <- 
  fung_mean_coreness %>% 
  full_join(bact_mean_coreness) %>% 
  full_join(full_mean_coreness)

### global efficency ####
fung_global_effic <- 
  data.frame(taxa = "Fungi",
             global_effic = I(
               complete_netstats$fungi %>% 
                 map('global_effic'))
  ) %>% mutate(sample_id = row.names(.))
bact_global_effic <- 
  data.frame(taxa = "Bacteria",
             global_effic = I(
               complete_netstats$bacteria %>% 
                 map('global_effic'))
  ) %>% mutate(sample_id = row.names(.))
full_global_effic <- 
  data.frame(taxa = "Both",
             global_effic = I(
               complete_netstats$full %>% 
                 map('global_effic'))
  ) %>% mutate(sample_id = row.names(.))
global_effic_df <- 
  fung_global_effic %>% 
  full_join(bact_global_effic) %>% 
  full_join(full_global_effic)

### clustering coeficient ####
fung_clustering_coeficent <- 
  data.frame(taxa = "Fungi",
             clustering_coeficent = I(
               complete_netstats$fungi %>% 
                 map('clustering_coeficient'))
  ) %>% mutate(sample_id = row.names(.))
bact_clustering_coeficent <- 
  data.frame(taxa = "Bacteria",
             clustering_coeficent = I(
               complete_netstats$bacteria %>% 
                 map('clustering_coeficient'))
  ) %>% mutate(sample_id = row.names(.))
full_clustering_coeficent <- 
  data.frame(taxa = "Both",
             clustering_coeficent = I(
               complete_netstats$full %>% 
                 map('clustering_coeficient'))
  ) %>% mutate(sample_id = row.names(.))
clustering_coeficent_df <- 
  fung_clustering_coeficent %>% 
  full_join(bact_clustering_coeficent) %>% 
  full_join(full_clustering_coeficent)
### Combine all attributes ####
graph_atributes_df <- 
  clique_num_df %>% 
  full_join(clustering_coeficent_df) %>% 
  full_join(deg_dist_df) %>% 
  full_join(global_effic_df) %>% 
  full_join(max_deg_df) %>% 
  full_join(mean_betweenness_df) %>% 
  full_join(mean_closeness_df) %>% 
  full_join(mean_coreness_df) %>% 
  full_join(mean_dist_df) %>% 
  full_join(n_vertices_df) %>% 
  rename('deg_distribution' = 'value')

# unlist list-columns
graph_atributes_df <- 
  graph_atributes_df %>% 
  mutate(across(where(is.list),unlist))
  
## Export attributes ####
graph_atributes_df %>% 
  saveRDS("./output/complete_network_attributes.RDS")

# PLOT NETWORK ATTRIBUTES ####

## Richness vs network ####
# get grouped richness estimates
data.frame(fung@sam_data$sample_id,fung@sam_data$treatment)
richness <- 
  fung %>% 
  merge_samples('sample_id') %>% 
  estimate_richness(measures = "Observed") %>% 
  mutate(taxa="Fungi", treatment = case_when(grepl(pattern="^P1",row.names(.)) ~ "1-2 years",
                                             grepl(pattern="^P2",row.names(.)) ~ "2-3 years",
                                             grepl(pattern="^P3",row.names(.)) ~ "3-4 years",
                                             grepl(pattern="^P4",row.names(.)) ~ "0-1 years",
                                             grepl(pattern="^P5",row.names(.)) ~ "oil palm"),
         sampleID = row.names(.)) %>% 
  full_join(
    bact %>% 
      merge_samples('sample_id') %>% 
      estimate_richness(measures = "Observed") %>% 
      mutate(taxa="Bacteria", treatment = case_when(grepl(pattern="^P1",row.names(.)) ~ "1-2 years",
                                                    grepl(pattern="^P2",row.names(.)) ~ "2-3 years",
                                                    grepl(pattern="^P3",row.names(.)) ~ "3-4 years",
                                                    grepl(pattern="^P4",row.names(.)) ~ "0-1 years",
                                                    grepl(pattern="^P5",row.names(.)) ~ "oil palm"),
             sampleID = row.names(.))
    ) %>% 
  full_join(
    full %>% 
      merge_samples('sample_id') %>% 
      estimate_richness(measures = "Observed") %>% 
      mutate(taxa="Both", treatment = case_when(grepl(pattern="^P1",row.names(.)) ~ "1-2 years",
                                                grepl(pattern="^P2",row.names(.)) ~ "2-3 years",
                                                grepl(pattern="^P3",row.names(.)) ~ "3-4 years",
                                                grepl(pattern="^P4",row.names(.)) ~ "0-1 years",
                                                grepl(pattern="^P5",row.names(.)) ~ "oil palm"),
             sampleID = row.names(.))
  ) 
richness_vs_atributes_df <- 
graph_atributes_df %>%
  rename("sampleID"="sample_id") %>% # rename temporarily for easier pivot
  dplyr::select(-deg_distribution) %>% 
  pivot_longer(contains("_"),names_to = "network_attribute") %>%
  full_join(richness) %>% 
  mutate(treatment = factor(treatment,levels=levels(fung@sam_data$treatment)),
         taxa=factor(taxa,levels=c("Bacteria","Fungi","Both"))) %>% 
  rename("sample_id" = "sampleID") # change name back

p1 <-   
richness_vs_atributes_df %>% 
  dplyr::filter(taxa=="Bacteria" & network_attribute != "clique_num") %>% 
  ggplot(aes(x=value,y=Observed,color=treatment)) +
  geom_point() +
  # geom_smooth(method='lm',se=FALSE) +
  facet_grid(cols=vars(network_attribute),scales = 'free') +
  scale_color_manual(values = trt.pal) +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5,size=10),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust=.5),
        legend.position = 'none') +
  labs(color="Taxa",y="",x="",title="Bacteria")

p2 <-   
  richness_vs_atributes_df %>% 
  dplyr::filter(taxa=="Fungi" & network_attribute != "clique_num") %>% 
  ggplot(aes(x=value,y=Observed,color=treatment)) +
  geom_point() +
  # geom_smooth(method='lm',se=FALSE) +
  facet_grid(cols=vars(network_attribute),scales = 'free') +
  scale_color_manual(values = trt.pal) +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5,size=10),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust=.5)) +
  labs(color="Treatment",y="Species richness",x="",title="Fungi")

p3 <-   
  richness_vs_atributes_df %>% 
  dplyr::filter(taxa=="Both" & network_attribute != "clique_num") %>% 
  ggplot(aes(x=value,y=Observed,color=treatment)) +
  geom_point() +
  # geom_smooth(method='lm',se=FALSE) +
  facet_grid(cols=vars(network_attribute),scales = 'free') +
  scale_color_manual(values = trt.pal) +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5,size=10),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust=.5),
        legend.position = 'none') +
  labs(color="Taxa",y="",x="",title = "Both") 

p1 / p2 / p3
### export plot ####
ggsave("./output/figs/Species_richness_v_network_attributes.png",dpi=400,height = 8,width = 16)

## treatment vs network ####
p4 <- richness_vs_atributes_df %>% 
  dplyr::filter(taxa=="Bacteria" & network_attribute != "clique_num") %>% 
  ggplot(aes(x=treatment,y=value,fill=treatment)) +
  geom_boxplot() +
  facet_wrap(~network_attribute,scales = 'free',nrow = 1) +
  labs(fill="Treatment",x="",y="",title="Bacteria") +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5,size=10),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust=.5)) +
  scale_fill_manual(values = trt.pal)
p5 <- richness_vs_atributes_df %>% 
  dplyr::filter(taxa=="Fungi" & network_attribute != "clique_num") %>% 
  ggplot(aes(x=treatment,y=value,fill=treatment)) +
  geom_boxplot() +
  facet_wrap(~network_attribute,scales = 'free',nrow = 1) +
  labs(fill="Treatment",x="",y="",title="Fungi") +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5,size=10),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust=.5)) +
  scale_fill_manual(values = trt.pal)
p6 <- richness_vs_atributes_df %>% 
  dplyr::filter(taxa=="Both" & network_attribute != "clique_num") %>% 
  ggplot(aes(x=treatment,y=value,fill=treatment)) +
  geom_boxplot() +
  facet_wrap(~network_attribute,scales = 'free',nrow = 1) +
  labs(fill="Treatment",x="",y="",title="Both") +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5,size=10),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust=.5)) +
  scale_fill_manual(values = trt.pal)

# combine plots
p4 / p5 / p6 + patchwork::plot_layout(guides = 'collect')

### export plot ####
ggsave("./output/figs/Treatment_v_Network_attributes.png",dpi=400,height = 8,width = 16)

# MODEL NETWORK ATTRIBUTES ####
glimpse(richness_vs_atributes_df)

## re-scale predictor values ####
scaled_richness_vs_atributes_df <- 
richness_vs_atributes_df %>% 
  rename('sampleID'='sample_id') %>% 
  pivot_wider(names_from = network_attribute,values_from = value,
              values_fn = mean) %>% 
  mutate(across(contains("_"),scale)) %>% 
  pivot_longer(contains("_"),names_to = 'network_attribute',values_to = 'value') %>% 
  mutate(value=as.numeric(value),
         Observed=as.numeric(Observed))

## richness ####
network_mod_fung <- 
  lmerTest::lmer(data=scaled_richness_vs_atributes_df %>% dplyr::filter(taxa=="Fungi"),
                 formula = Observed ~ value + value:network_attribute + (1| treatment))
network_mod_bact <- 
  lmerTest::lmer(data=scaled_richness_vs_atributes_df %>% dplyr::filter(taxa=="Bacteria"),
                 formula = Observed ~ value + value:network_attribute + (1| treatment))
network_mod_both <- 
  lmerTest::lmer(data=scaled_richness_vs_atributes_df %>% dplyr::filter(taxa=="Both"),
                 formula = Observed ~ value + value:network_attribute + (1| treatment))

## export model results ####

lmerTest:::get_coefmat(network_mod_bact) %>% 
  as.data.frame() %>% 
  mutate(term=row.names(.) %>% str_remove('network_attribute'),
         taxa="Bacteria") %>% 
  full_join(
    lmerTest:::get_coefmat(network_mod_fung) %>% 
      as.data.frame() %>% 
      mutate(term=row.names(.) %>% str_remove('network_attribute'),
             taxa="Fungi")
  ) %>% 
  full_join(
    lmerTest:::get_coefmat(network_mod_both) %>% 
      as.data.frame() %>% 
      mutate(term=row.names(.) %>% str_remove('network_attribute'),
             taxa="Both")
  ) %>% 
  dplyr::select(taxa,term,Estimate,`Std. Error`,df,`t value`,`Pr(>|t|)`) %>% 
  write_csv("./output/Network_vs_richness_model_table.csv")

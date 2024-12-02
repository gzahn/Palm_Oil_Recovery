# COMMUNITY ASSEMBLY OVER TIME

# SETUP ####

## packages ####
library(tidyverse)
library(phyloseq)
library(janitor)
library(patchwork)
library(broom)
library(vegan)
library(ShortRead)

## functions ####
source("./R/functions.R")

## options ####
theme_set(theme_bw() +
            theme(axis.title = element_text(face='bold',size=18),
                  axis.text = element_text(face='bold',size=14),
                  legend.title = element_text(face='bold',size=18),
                  legend.text = element_text(face='bold',size=14),
                  strip.text = element_text(face='bold',size=18),
                  plot.title = element_text(face='bold',size=18))
)

set.seed(666)


## load physeq objects ####
fung <- readRDS("./output/ITS_Physeq_cleaned_w_tree.RDS")
bact <- readRDS("./output/16S_Physeq_cleaned_w_tree.RDS")
# add refseqs so we can rename the taxa
fung@refseq <- taxa_names(fung) %>% DNAStringSet()
bact@refseq <- taxa_names(bact) %>% DNAStringSet()
taxa_names(fung) <- paste0("ASV_",1:ntaxa(fung))
taxa_names(bact) <- paste0("ASV_",1:ntaxa(bact))

# add numeric "time since palm oil" variable
fung@sam_data$time <- 
  fung@sam_data$treatment %>% 
  as.character() %>% 
  str_sub(1,1) %>% 
  str_replace("3",'4') %>% 
  str_replace("2",'3') %>% 
  str_replace("1",'2') %>% 
  str_replace("0",'1') %>% 
  str_replace("o",'0') %>% 
  as.numeric()

bact@sam_data$time <- 
  bact@sam_data$treatment %>% 
  as.character() %>% 
  str_sub(1,1) %>% 
  str_replace("3",'4') %>% 
  str_replace("2",'3') %>% 
  str_replace("1",'2') %>% 
  str_replace("0",'1') %>% 
  str_replace("o",'0') %>% 
  as.numeric()

# View tree
bact_order.glom <- 
  bact %>% 
  tax_glom("Order")

plot_tree(bact_order.glom, color="Phylum")
# need a tree of ALL bacterial orders:
# highlight the diversity (at order level) present at each year, from left to right

# agglomerate taxa at species level (not ASV level for this analysis)
fung <- 
  fung %>% 
  tax_glom("Species",NArm = FALSE,bad_empty = c(NA, "", " ", "\t","sp."))
bact <- 
  bact %>% 
  tax_glom("Species",NArm = FALSE,bad_empty = c(NA, "", " ", "\t","sp."))

# find basic pool of taxa present in palm oil plots
fung_t0 <- 
  fung %>% 
  subset_samples(time == 0)
fung_t0 <- 
  fung_t0 %>% 
  subset_taxa(taxa_sums(fung_t0) > 0)
fung_t1 <- 
  fung %>% 
  subset_samples(time == 1)
fung_t1 <- 
  fung_t1 %>% 
  subset_taxa(taxa_sums(fung_t1) > 0)
fung_t2 <- 
  fung %>% 
  subset_samples(time == 2)
fung_t2 <- 
  fung_t2 %>% 
  subset_taxa(taxa_sums(fung_t2) > 0)
fung_t3 <- 
  fung %>% 
  subset_samples(time == 3)
fung_t3 <- 
  fung_t3 %>% 
  subset_taxa(taxa_sums(fung_t3) > 0)
fung_t4 <- 
  fung %>% 
  subset_samples(time == 4)
fung_t4 <- 
  fung_t4 %>% 
  subset_taxa(taxa_sums(fung_t4) > 0)

bact_t0 <- 
  bact %>% 
  subset_samples(time == 0)
bact_t0 <- 
  bact_t0 %>% 
  subset_taxa(taxa_sums(bact_t0) > 0)
bact_t1 <- 
  bact %>% 
  subset_samples(time == 1)
bact_t1 <- 
  bact_t1 %>% 
  subset_taxa(taxa_sums(bact_t1) > 0)
bact_t2 <- 
  bact %>% 
  subset_samples(time == 2)
bact_t2 <- 
  bact_t2 %>% 
  subset_taxa(taxa_sums(bact_t2) > 0)
bact_t3 <- 
  bact %>% 
  subset_samples(time == 3)
bact_t3 <- 
  bact_t3 %>% 
  subset_taxa(taxa_sums(bact_t3) > 0)
bact_t4 <- 
  bact %>% 
  subset_samples(time == 4)
bact_t4 <- 
  bact_t4 %>% 
  subset_taxa(taxa_sums(bact_t4) > 0)

# which taxa in each year after palm oil are members of the palm oil community?
# which taxa are "new"?


# which taxa are "continually new" .. as in, not found in any previous year ?
fung_t1 <- 
  fung_t1 %>% 
  subset_taxa(
    taxa_names(fung_t1) %ni% taxa_names(fung_t0)
)
fung_t2 <- 
  fung_t2 %>% 
  subset_taxa(
    taxa_names(fung_t2) %ni% c(taxa_names(fung_t0),taxa_names(fung_t1))
  )
fung_t3 <- 
  fung_t3 %>% 
  subset_taxa(
    taxa_names(fung_t3) %ni% c(taxa_names(fung_t0),taxa_names(fung_t1),taxa_names(fung_t2))
  )
fung_t4 <- 
  fung_t4 %>% 
  subset_taxa(
    taxa_names(fung_t4) %ni% c(taxa_names(fung_t0),taxa_names(fung_t1),taxa_names(fung_t2),taxa_names(fung_t3))
  )

bact_t1 <- 
  bact_t1 %>% 
  subset_taxa(
    taxa_names(bact_t1) %ni% taxa_names(bact_t0)
  )
bact_t2 <- 
  bact_t2 %>% 
  subset_taxa(
    taxa_names(bact_t2) %ni% c(taxa_names(bact_t0),taxa_names(bact_t1))
  )
bact_t3 <- 
  bact_t3 %>% 
  subset_taxa(
    taxa_names(bact_t3) %ni% c(taxa_names(bact_t0),taxa_names(bact_t1),taxa_names(bact_t2))
  )
bact_t4 <- 
  bact_t4 %>% 
  subset_taxa(
    taxa_names(bact_t4) %ni% c(taxa_names(bact_t0),taxa_names(bact_t1),taxa_names(bact_t2),taxa_names(bact_t3))
  )

# remove phylogenetic trees to recombine the objects
fung_t0@phy_tree <- NULL
fung_t1@phy_tree <- NULL
fung_t2@phy_tree <- NULL
fung_t3@phy_tree <- NULL
fung_t4@phy_tree <- NULL
bact_t0@phy_tree <- NULL
bact_t1@phy_tree <- NULL
bact_t2@phy_tree <- NULL
bact_t3@phy_tree <- NULL
bact_t4@phy_tree <- NULL

# how to plot this? ....
#     NEW TAXA THAT APPEAR OVER TIME?

# combine the objects
fung_timeseries <- merge_phyloseq(fung_t0,fung_t1,fung_t2,fung_t3,fung_t4)
bact_timeseries <- merge_phyloseq(bact_t0,bact_t1,bact_t2,bact_t3,bact_t4)
# remove empty samples
fung_timeseries <- 
  fung_timeseries %>% 
  subset_samples(sample_sums(fung_timeseries) > 0)
bact_timeseries <- 
  bact_timeseries %>% 
  subset_samples(sample_sums(bact_timeseries) > 0)



fung_timeseries %>% 
  merge_samples('plot')
fung_timeseries@sam_data$time <- factor(fung_timeseries@sam_data$time)
fung_timeseries %>% 
  # subset_samples(time != 0) %>%
  
  merge_samples('time') %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Order")



fung_timeseries %>% 
  transform_sample_counts(pa) %>% 
  plot_bar2(x='time',fill="Phylum") +
  labs(y="Richness",x="Years since palm oil")

# count number of unique taxa that are "arriving" over time in each sample (by year)


fung_timeseries %>% 
  transform_sample_counts(pa) %>% 
  psmelt() %>% 
  group_by(time,Sample) %>% 
  summarize(ntaxa=sum(Abundance)) %>% 
  mutate(group="Fungi") %>% 
full_join(
  bact_timeseries %>% 
    transform_sample_counts(pa) %>% 
    psmelt() %>% 
    group_by(time,Sample) %>% 
    summarize(ntaxa=sum(Abundance)) %>% 
    mutate(group="Bacteria")
  ) %>% 
  write_csv("./output/taxa_counts_by_sample_over_time.csv")

# total number of new taxa in each year (not by site)
data.frame(
  year=rep(0:4,2),
  taxon=rep(c("Fungi","Bacteria"),each=5),
  n_new_taxa=c(ntaxa(fung_t0),ntaxa(fung_t1),ntaxa(fung_t2),ntaxa(fung_t3),ntaxa(fung_t4),
               ntaxa(bact_t0),ntaxa(bact_t1),ntaxa(bact_t2),ntaxa(bact_t3),ntaxa(bact_t4)),
  resolution="species"
) %>% 
  write_csv("./output/total_new_taxa_over_time.csv")

# Just exactly WHO is arriving over time?

# build data frames of relabund info
fung_ts_melt <- 
fung_timeseries %>% 
  transform_sample_counts(ra) %>% 
  psmelt()
bact_ts_melt <- 
  bact_timeseries %>% 
  transform_sample_counts(ra) %>% 
  psmelt()
# add full taxonomy column
fung_ts_melt <- 
  fung_ts_melt %>% 
  mutate(taxonomy=paste(Phylum,Class,Order,Family,Genus,Species))
bact_ts_melt <- 
  bact_ts_melt %>% 
  mutate(taxonomy=paste(Phylum,Class,Order,Family,Genus,Species))


# find taxa that have >= 5% relative abundance
# these are "new" taxa that are present in significant numbers
fung_new_arrivals <- 
fung_ts_melt %>% 
  group_by(time,OTU) %>% 
  reframe(pH=mean(pH),
          c_perc_ww=mean(c_perc_ww),
          h_perc_ww=mean(h_perc_ww),
          n_perc_ww=mean(n_perc_ww),
          p_ppm=mean(p_ppm),
          mean_relabund=mean(Abundance),
          Taxonomy=unique(taxonomy) %>% str_remove(" NA") %>% str_remove("\\(Fungi\\)")) %>% 
  dplyr::filter(mean_relabund >= .05 & time != 0) %>% 
  mutate(Taxonomy=factor(Taxonomy,levels=Taxonomy))


bact_new_arrivals <- 
bact_ts_melt %>% 
  group_by(time,OTU) %>% 
  reframe(pH=mean(pH),
          c_perc_ww=mean(c_perc_ww),
          h_perc_ww=mean(h_perc_ww),
          n_perc_ww=mean(n_perc_ww),
          p_ppm=mean(p_ppm),
          mean_relabund=mean(Abundance),
          Taxonomy=unique(taxonomy) %>% str_remove(" NA")) %>% 
  dplyr::filter(mean_relabund >= .05 & time != 0) %>% 
  mutate(Taxonomy=factor(Taxonomy,levels=Taxonomy))

fung_new_arrivals %>% 
  mutate(group="Fungi",
         time=as.numeric(time)) %>% 
  full_join(
    bact_new_arrivals %>% 
      mutate(group="Bacteria")
  ) %>% 
  write_csv("./output/new_arrival_taxa_over_time.csv")


p1 <- 
fung_new_arrivals %>% 
  ggplot(aes(x=mean_relabund,y=fct_rev(Taxonomy),
             color=as.character(time),
             fill=factor(time))) +
  geom_col(color='white') +
  scale_fill_manual(values=trt.pal[2:5]) +
  theme(axis.text.y = element_text(face='bold.italic',size=10),
        legend.position = "top") +
  labs(x="Mean relative abundance",y="Fungi",fill="Years since restoration")

p2 <- 
bact_new_arrivals %>% 
  ggplot(aes(x=mean_relabund,y=fct_rev(Taxonomy),
             color=as.character(time),
             fill=factor(time))) +
  geom_col(color='white') +
  scale_fill_manual(values=trt.pal[2:5]) +
  theme(axis.text.y = element_text(face='bold.italic',size=10),
        legend.position = 'none') +
  labs(x="Mean relative abundance",y="Bacteria",fill="Years since restoration")

p1 / p2
ggsave("./output/figs/Newly_arriving_taxa_over_time_bact_and_fung.png",height = 12, width = 16,dpi=400)


# Export data for Carl ####

# reload raw data
fung <- readRDS("./output/ITS_Physeq_cleaned_w_tree.RDS")
bact <- readRDS("./output/16S_Physeq_cleaned_w_tree.RDS")

# agglomerate taxa at order level, merge by treatment year, and find relabund
ord_bact <- 
  bact %>% 
  tax_glom("Order",NArm = FALSE) %>% 
  merge_samples("treatment") %>% 
  transform_sample_counts(ra)

ord_fung <- 
  fung %>% 
  tax_glom("Order",NArm = FALSE) %>% 
  merge_samples("treatment") %>% 
  transform_sample_counts(ra)



# melt objects to data frames
ord_bact_df <- 
  ord_bact %>% 
  psmelt() %>% 
  mutate(time = Sample %>% 
           as.character() %>% 
           str_sub(1,1) %>% 
           str_replace("3",'4') %>% 
           str_replace("2",'3') %>% 
           str_replace("1",'2') %>% 
           str_replace("0",'1') %>% 
           str_replace("o",'0') %>% 
           as.numeric())
ord_fung_df <- 
  ord_fung %>% 
  psmelt() %>% 
  mutate(time = Sample %>% 
           as.character() %>% 
           str_sub(1,1) %>% 
           str_replace("3",'4') %>% 
           str_replace("2",'3') %>% 
           str_replace("1",'2') %>% 
           str_replace("0",'1') %>% 
           str_replace("o",'0') %>% 
           as.numeric())

# remove rows where relabund == 0
# need to find number of unique otus for each subset (time,order)
b_melt <- 
  bact %>% 
  psmelt()
f_melt <- 
  fung %>% 
  psmelt()

bact_otu_counts <- 
  b_melt %>% 
  mutate(time = treatment %>% 
           as.character() %>% 
           str_sub(1,1) %>% 
           str_replace("3",'4') %>% 
           str_replace("2",'3') %>% 
           str_replace("1",'2') %>% 
           str_replace("0",'1') %>% 
           str_replace("o",'0') %>% 
           as.numeric()) %>% 
  dplyr::filter(Abundance > 0) %>% 
  group_by(Order,time) %>% 
  summarize(N_ASVs = n())

fung_otu_counts <- 
  f_melt %>% 
  mutate(time = treatment %>% 
           as.character() %>% 
           str_sub(1,1) %>% 
           str_replace("3",'4') %>% 
           str_replace("2",'3') %>% 
           str_replace("1",'2') %>% 
           str_replace("0",'1') %>% 
           str_replace("o",'0') %>% 
           as.numeric()) %>% 
  dplyr::filter(Abundance > 0) %>% 
  group_by(Order,time) %>% 
  summarize(N_ASVs = n())

ord_bact_df %>% 
  dplyr::rename("ASV" = "OTU",
                "relative_abundance" = "Abundance") %>% 
  dplyr::select(Kingdom,Phylum,Class,Order,time,relative_abundance,ASV) %>% 
  dplyr::filter(!is.na(Order) & relative_abundance > 0 & Order != "unclassified") %>% 
  dplyr::arrange(Order,time) %>%
  left_join(bact_otu_counts) %>% 
  write_csv("./output/Bacterial_order_presence_over_time.csv")


ord_fung_df %>% 
  dplyr::rename("ASV" = "OTU",
                "relative_abundance" = "Abundance") %>% 
  dplyr::select(Kingdom,Phylum,Class,Order,time,relative_abundance,ASV) %>% 
  dplyr::filter(!is.na(Order) & relative_abundance > 0 & Order != "unclassified") %>%
  dplyr::arrange(Order,time) %>% 
  left_join(fung_otu_counts) %>% 
  write_csv("./output/Fungal_order_presence_over_time.csv")

x <- ord_fung_df %>% 
  dplyr::rename("ASV" = "OTU",
                "relative_abundance" = "Abundance") %>% 
  dplyr::select(Kingdom,Phylum,Class,Order,time,relative_abundance,ASV) %>% 
  dplyr::filter(!is.na(Order) & relative_abundance > 0 & Order != "unclassified") %>%
  dplyr::arrange(Order,time) %>% 
  left_join(fung_otu_counts)

x %>% 
  group_by(time) %>% 
  summarize(total_abund = sum(relative_abundance))

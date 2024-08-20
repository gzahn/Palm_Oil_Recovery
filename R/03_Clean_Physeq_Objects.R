# CLEAN PHYSEQ OBJECTS

# SETUP ####
set.seed(666)

## packages ####
library(tidyverse)
library(dada2)
library(phyloseq)
library(decontam)
library(janitor)
library(patchwork)

## functions ####
source("./R/functions.R")

## load metadata ####
meta <- read_csv("./data/Metadata.csv") %>% 
  janitor::clean_names()

## load physeq objects ####
ps_its <- readRDS("./output/ITS_Physeq_not-cleaned.RDS")
ps_16s <- readRDS("./output/16S_Physeq_not-cleaned.RDS")


# EXAMINE POSITIVE CONTROLS ####

# subset to pos controls
pos_16s <- ps_16s %>% 
  subset_samples(control == "Pos")
pos_16s <- pos_16s %>% 
  subset_taxa(taxa_sums(pos_16s) > 0)
pos_its <- ps_its %>% 
  subset_samples(control == "Pos")
pos_its <- pos_its %>% 
  subset_taxa(taxa_sums(pos_its) > 0)

# barplots
pos_16s %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Genus")
pos_its %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Genus")

# Zymo mock community members
zymo <- c("Pseudomonas aeruginosa","Escherichia coli","Salmonella enterica","Lactobacillus fermentum",
  "Enterococcus faecalis","Staphylococcus aureus","Listeria monocytogenes","Bacillus subtilis",
  "Saccharomyces cerevisiae","Cryptococcus neoformans")
zymo_genera <- str_split(zymo," ") %>% map_chr(1)


# rename "Species"
pos_16s@tax_table[,7] <- 
paste(pos_16s@tax_table[,6],pos_16s@tax_table[,7]) %>% 
  str_replace_all(" NA"," sp.")
pos_its@tax_table[,7] <- 
  paste(pos_its@tax_table[,6],pos_its@tax_table[,7]) %>% 
  str_replace_all(" NA"," sp.")

pos_plot1 <- 
pos_16s %>% 
  subset_taxa(pos_16s %>% 
                transform_sample_counts(ra) %>% 
                otu_table() %>% colMeans() > .01) %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Species") + scale_fill_viridis_d(option='turbo')
pos_plot2 <- 
pos_its %>% 
  subset_taxa(pos_its %>% 
                transform_sample_counts(ra) %>% 
                otu_table() %>% colMeans() > .01) %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Species") + scale_fill_viridis_d(option='turbo')

pos_plot1 | pos_plot2
ggsave("./output/figs/positive_control_barplots.png",width = 8,height = 6)

# Looks good...

# REMOVE POSITIVE CONTROLS ####
ps_16s <- ps_16s %>% 
  subset_samples(is.na(control))
ps_16s <- ps_16s %>% 
  subset_taxa(taxa_sums(ps_16s) > 1)
ps_its <- ps_its %>% 
  subset_samples(is.na(control))
ps_its <- ps_its %>% 
  subset_taxa(taxa_sums(ps_its) > 1)

# REMOVE NON_TARGET TAXA ####
ps_16s <- 
  ps_16s %>% 
  subset_taxa(!is.na(Phylum) & Family != "Mitochondria" & Order != "Chloroplast")
ps_its <- 
  ps_its %>% 
  subset_taxa(Kingdom == "k__Fungi" & !is.na(Phylum))



# CLEAN TAXONOMIC ASSIGNMENTS ####
ps_its <- 
  ps_its %>% 
  clean_ps_taxonomy()
ps_its@tax_table[,7] <- 
  ifelse(as.character(ps_its@tax_table[,7]) == "unclassified", NA, as.character(ps_its@tax_table[,7])) %>% 
  ifelse(is.na(.),"sp.",.)
ps_16s@tax_table[,7] <- 
  ifelse(as.character(ps_16s@tax_table[,7]) == "unclassified", NA, as.character(ps_16s@tax_table[,7])) %>% 
  ifelse(is.na(.),"sp.",.)


# REMOVE LOW-COUNT TAXA ####
ps_16s <- 
  ps_16s %>% 
  subset_taxa(taxa_sums(ps_16s) > 99)

# CLEAN UP METADATA ####
cols_to_keep <- c("sample_name","sample_id","p_h","c_perc_ww","h_perc_ww","n_perc_ww","s_perc_ww","p_ppm",
                  "plot","treatment","latitude","longitude","planting_regime","region")

ps_16s@sam_data <- ps_16s@sam_data[,cols_to_keep]
ps_its@sam_data <- ps_its@sam_data[,cols_to_keep]

glimpse(as(ps_16s@sam_data,'data.frame'))

# coerce column types
  # empty column
ps_16s@sam_data$s_perc_ww <- NULL
ps_its@sam_data$s_perc_ww <- NULL
  # treatment order
ps_16s@sam_data$treatment <- 
  ps_16s@sam_data$treatment %>% 
  factor(levels=c("oil palm","0-1 years","1-2 years","2-3 years","3-4 years"))
ps_its@sam_data$treatment <- 
  ps_its@sam_data$treatment %>% 
  factor(levels=c("oil palm","0-1 years","1-2 years","2-3 years","3-4 years"))
  # 'plot' to factor
ps_16s@sam_data$plot <- 
  ps_16s@sam_data$plot %>% factor()
ps_its@sam_data$plot <- 
  ps_its@sam_data$plot %>% factor()
  # rename pH
names(ps_16s@sam_data)[names(ps_16s@sam_data) == "p_h"] <- "pH"
names(ps_its@sam_data)[names(ps_its@sam_data) == "p_h"] <- "pH"

# EXPORT CLEAN PHYSEQ OBJECTS ####
saveRDS(ps_16s,"./output/16S_Physeq_cleaned.RDS")
saveRDS(ps_its,"./output/ITS_Physeq_cleaned.RDS")

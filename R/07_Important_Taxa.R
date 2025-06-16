# IMPORANT TAXA AND DIFFERENTIAL ABUNDANCE

# SETUP ####

## packages ####
library(tidyverse)
library(phyloseq)
library(janitor)
library(patchwork)
library(broom)
library(vegan)
library(corncob);packageVersion('corncob')
library(ShortRead)
library(vip)
library(ranger)

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
fung <- readRDS("./output/ITS_Physeq_cleaned_w_tree.RDS") %>% tax_glom("Species")
bact <- readRDS("./output/16S_Physeq_cleaned_w_tree.RDS") %>% tax_glom("Species")

# add refseqs so we can rename the taxa
fung@refseq <- taxa_names(fung) %>% DNAStringSet()
bact@refseq <- taxa_names(bact) %>% DNAStringSet()
taxa_names(fung) <- paste0("ASV_",1:ntaxa(fung))
taxa_names(bact) <- paste0("ASV_",1:ntaxa(bact))




# DIFFERENTIAL TEST ####

## fungi ####
fung_da_analysis <- differentialTest(formula = ~ treatment, #abundance
                                     phi.formula = ~ 1, #dispersion
                                     formula_null = ~ 1, #mean
                                     phi.formula_null = ~ 1,
                                     test = "Wald", boot = FALSE,
                                     data = fung,
                                     fdr_cutoff = 0.05,
                                     full_output = TRUE)
plot(fung_da_analysis)
saveRDS(fung_da_analysis,"./output/fungi_diff_analysis.RDS")
fung_da_analysis <- readRDS("./output/fungi_diff_analysis.RDS")
# name model objects
names(fung_da_analysis$full_output) <- taxa_names(fung)
# extract models for significant taxa
fung_sig_models <- fung_da_analysis$full_output[fung_da_analysis$significant_taxa]


## bacteria ####
bact_da_analysis <- differentialTest(formula = ~ treatment, #abundance
                                     phi.formula = ~ 1, #dispersion
                                     formula_null = ~ 1, #mean
                                     phi.formula_null = ~ 1,
                                     test = "Wald", boot = FALSE,
                                     data = bact,
                                     fdr_cutoff = 0.05,
                                     full_output = TRUE)
plot(bact_da_analysis)
saveRDS(bact_da_analysis,"./output/bacteria_diff_analysis.RDS")
bact_da_analysis <- readRDS("./output/bacteria_diff_analysis.RDS")
# name model objects
names(bact_da_analysis$full_output) <- taxa_names(bact)
# extract models for significant taxa
bact_sig_models <- bact_da_analysis$full_output[bact_da_analysis$significant_taxa]

bact_sig_models[[3]]
# REFINE SIG TAXA ####

## increased taxa ####

# fungi
# find models where estimate is always greater than intercept (palm oil site)
fung_greater <- logical()
for(mod in names(fung_sig_models)){
  x <- fung_sig_models[[mod]][['b.mu']]
  fung_greater[mod] <- is_greater_than_first(x)
}
# subset to only true values
fung_greater <- fung_greater[fung_greater == TRUE]


# bacteria
# find models where estimate is always greater than intercept (palm oil site)
bact_greater <- logical()
for(mod in names(bact_sig_models)){
  x <- bact_sig_models[[mod]][['b.mu']]
  bact_greater[mod] <- is_increasing(x)
}
# subset to only true values
bact_greater <- bact_greater[bact_greater == TRUE]
# can subset to just those taxa to keep the story simple (not simplistic)
names(fung_greater)



## randforest mods ####

# fungi
fungi_rf_mod <- 
data.frame(treatment = fung@sam_data$treatment) %>% 
  bind_cols(fung %>% 
              transform_sample_counts(ra) %>% 
              otu_table() %>% 
              as('matrix')) %>% 
  ranger(formula = treatment ~ .,data=., importance = 'permutation')
# most important taxa (importance > 3rd quantile of importance values)
fungi_top_taxa <- vi(fungi_rf_mod) %>% 
  dplyr::filter(Importance > quantile(Importance,.75)) %>% 
  pluck("Variable")

# find the significant taxa that increased over time, subset to them, melt
# merge by treatment
fung_sig_melt <- fung 
  # fung %>% 
  # merge_samples('treatment')
fung_sig_melt@sam_data$treatment <- factor(c("oil palm","0-1 years","1-2 years","2-3 years","3-4 years"),
                                           levels=c("oil palm","0-1 years","1-2 years","2-3 years","3-4 years"))
fung_sig_melt <- 
  fung_sig_melt %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(taxa_names(fung) %in% names(fung_greater) &
                taxa_names(fung) %in% fungi_top_taxa) %>% 
  psmelt() %>% 
  mutate(taxonomy=paste(Phylum,Order,Family,Genus,Species,sep=" ") %>% str_remove_all("NA "))

# bacteria
bacteria_rf_mod <- 
  data.frame(treatment = bact@sam_data$treatment) %>% 
  bind_cols(bact %>% 
              transform_sample_counts(ra) %>% 
              otu_table() %>% 
              as('matrix')) %>% 
  ranger(formula = treatment ~ .,data=., importance = 'permutation')
# most important taxa (importance > 3rd quantile of importance values)
bacteria_top_taxa <- vi(bacteria_rf_mod) %>% 
  dplyr::filter(Importance > quantile(Importance,.75)) %>% 
  pluck("Variable")

# find the significant taxa that increased over time, subset to them, melt
# merge by treatment
bact_sig_melt <- bact 
# bact %>% 
# merge_samples('treatment')
bact_sig_melt@sam_data$treatment <- factor(c("oil palm","0-1 years","1-2 years","2-3 years","3-4 years"),
                                           levels=c("oil palm","0-1 years","1-2 years","2-3 years","3-4 years"))
bact_sig_melt <- 
  bact_sig_melt %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(taxa_names(bact) %in% names(bact_greater) &
                taxa_names(bact) %in% bacteria_top_taxa) %>% 
  psmelt() %>% 
  mutate(taxonomy=paste(Phylum,Order,Family,Genus,Species,sep=" ") %>% str_remove_all("NA "))

saveRDS(fung_sig_melt,"./output/fungi_significant_taxa.RDS")
saveRDS(bact_sig_melt,"./output/bacteria_significant_taxa.RDS")


names(bact_greater)
bact_sig_melt$OTU %in% names(bact_greater)
bact_sig_melt %>% 
  group_by(treatment,taxonomy) %>% 
  summarize(relabund = mean(Abundance,na.rm=TRUE))

  



# PLOTS ####

## Plot important taxa relative abundances by timepoint ####
p_fung_imp_taxa <- 
fung %>% 
  merge_samples("treatment",fun = 'sum') %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(taxa_names(fung) %in% unique(fung_sig_melt$OTU)) %>% 
  psmelt() %>% 
  mutate(Sample = factor(Sample,levels=levels(fung@sam_data$treatment))) %>% 
  ggplot(aes(x=Sample,y=Abundance,fill=Class)) +
  geom_col() +
  scale_fill_viridis_d(option = 'turbo') + theme(plot.caption = element_text(face='bold',size=12)) +
  labs(title="Fungi",y="Relative abundance",x="\nRewilding stage",caption = paste("N taxa = ",length(unique(fung_sig_melt$OTU))))

p_bact_imp_taxa <- 
bact %>% 
  merge_samples("treatment",fun = 'sum') %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(taxa_names(bact) %in% unique(bact_sig_melt$OTU)) %>% 
  psmelt() %>% 
  mutate(Sample = factor(Sample,levels=levels(bact@sam_data$treatment))) %>% 
  ggplot(aes(x=Sample,y=Abundance,fill=Class)) +
  geom_col() +
  scale_fill_viridis_d(option = 'turbo') + theme(plot.caption = element_text(face='bold',size=12)) +
  labs(title="Bacteria",y="Relative abundance",x="",caption = paste("N taxa = ",length(unique(bact_sig_melt$OTU))," "))

p_bact_imp_taxa / p_fung_imp_taxa
ggsave("./output/figs/important_taxa_by_time.png",dpi=400,height = 9,width = 12)

bact_sig_melt %>% 
  # dplyr::filter(treatment != "oil palm") %>% 
  ggplot(aes(x=treatment,y=Abundance,fill=Genus)) +
  geom_boxplot() +
  facet_wrap(~Genus,scales = 'free') +
  theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)) +
  scale_fill_viridis_d(option="turbo")
ggsave("./output/figs/important_bacterial_taxa.png",height = 16, width = 16, dpi=400)

fung_sig_melt %>% 
  # dplyr::filter(treatment != "oil palm") %>% 
  ggplot(aes(x=treatment,y=Abundance,fill=Order)) +
  geom_boxplot() +
  facet_wrap(~Order,scales = 'free') +
  theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)) +
  scale_fill_viridis_d(option="turbo")
ggsave("./output/figs/important_fungal_taxa.png",height = 16, width = 16, dpi=400)



# bact_sig_melt <- 
#   bact %>% 
#   transform_sample_counts(ra) %>% 
#   subset_taxa(taxa_names(bact) %in% bact_da_analysis$significant_taxa) %>% 
#   psmelt()
# 
# 
# 
# fung_da_analysis$significant_taxa
# fung_da_analysis$significant_taxa %>% otu_to_taxonomy(data=fung)

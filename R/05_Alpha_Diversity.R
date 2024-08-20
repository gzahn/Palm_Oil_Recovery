# ALPHA DIVERSITY 

# SETUP ####
set.seed(666)

## packages ####
library(tidyverse)
library(dada2)
library(phyloseq)
library(decontam)
library(janitor)
library(patchwork)
library(GGally)
library(lmerTest)
library(ggpubr)
library(broom)
library(easystats)
library(multcompView)

## functions ####
source("./R/functions.R")

## options ####
theme_set(theme_bw() +
            theme(axis.title = element_text(face='bold',size=18),
                  axis.text = element_text(face='bold',size=14),
                  legend.title = element_text(face='bold',size=18),
                  legend.text = element_text(face='bold',size=14),
                  strip.text = element_text(face='bold',size=18))
)

## load physeq objects ####
fung <- readRDS("./output/ITS_Physeq_cleaned_w_tree.RDS")
bact <- readRDS("./output/16S_Physeq_cleaned_w_tree.RDS")

## prep external metadata ####
meta_fung <- fung@sam_data %>% as('data.frame')
meta_bact <- bact@sam_data %>% as('data.frame')
full_meta <- full_join(meta_fung,meta_bact)
meta <- meta_fung %>% 
  dplyr::select(pH,c_perc_ww,h_perc_ww,n_perc_ww,p_ppm,treatment,plot,sample_id, latitude,longitude)
# remove "-" from treatment for multcomView
meta$treatment <- 
  meta$treatment %>% 
  as.character() %>% str_replace("-","_") %>% 
  factor(levels=c("oil palm","0_1 years","1_2 years","2_3 years","3_4 years"))

### calculate richness ####
# ASV richness
meta$fungal_asv_richness <- 
  fung %>% 
  estimate_richness(measures = 'Observed') %>% 
  pluck('Observed')
meta$bacterial_asv_richness <- 
  bact %>% 
  estimate_richness(measures = 'Observed') %>% 
  pluck('Observed')
# taxon (species) richness
meta$fungal_spp_richness <- 
  fung %>% 
  tax_glom("Species") %>% 
  estimate_richness(measures = 'Observed') %>% 
  pluck('Observed')
meta$bacterial_spp_richness <- 
  bact %>% 
  tax_glom("Species") %>% 
  estimate_richness(measures = 'Observed') %>% 
  pluck('Observed')


# BASIC READ STATS ####
asv_counts <- 
  data.frame(
    fungal_asvs = phyloseq::ntaxa(fung),
    fungal_spp = fung %>% 
      tax_glom("Species") %>% 
      phyloseq::ntaxa(),
    fungal_gen = fung %>% 
      tax_glom("Genus") %>% 
      phyloseq::ntaxa(),
    bacterial_asvs = phyloseq::ntaxa(bact),
    bacterial_spp = bact %>% 
      tax_glom("Species") %>% 
      phyloseq::ntaxa(),
    bacterial_gen = bact %>% 
      tax_glom("Genus") %>% 
      phyloseq::ntaxa()
  )
write_csv(asv_counts,"./output/ASV_and_spp_counts.csv")


# EXPLORE METADATA ####
full_meta %>% 
  dplyr::select(pH,c_perc_ww,h_perc_ww,n_perc_ww,p_ppm,treatment) %>% 
  ggpairs(mapping=ggplot2::aes(color=treatment)) +
  scale_color_manual(values = trt.pal) +
  scale_fill_manual(values = trt.pal)

table(meta$plot,meta$treatment)



# METADATA MODELS ####

## full model ####
meta %>% 
  pivot_longer(all_of(c("pH","c_perc_ww","h_perc_ww","n_perc_ww","p_ppm")),
               names_to = 'Attribute',values_to = "Value") %>% 
  mutate(Attribute = Attribute %>% 
           str_to_sentence() %>% 
           str_replace("_perc_ww"," (% ww)") %>% 
           str_replace("_ppm"," (ppm)") %>% 
           str_replace("Ph","pH")) %>% 
  ggplot(aes(x=treatment,y=Value,fill=treatment)) +
  geom_boxplot() +
  facet_wrap(~Attribute,scales = 'free') +
  scale_fill_manual(values = trt.pal) +
  labs(x="Treatment",fill="Treatment")

## Individual ANOVA models ####
# pH
mod_pH <- 
  aov(data = meta,
      formula = pH ~ treatment)
broom::tidy(mod_pH) %>% 
  write_csv("./output/pH_vs_treatment_aov-model.csv")
# C
mod_C <- 
  aov(data = meta,
      formula = c_perc_ww ~ treatment)
broom::tidy(mod_C) %>% 
  write_csv("./output/C_vs_treatment_aov-model.csv")
# H
mod_H <- 
  aov(data = meta,
      formula = h_perc_ww ~ treatment)
broom::tidy(mod_H) %>% 
  write_csv("./output/H_vs_treatment_aov-model.csv")
# N
mod_N <- 
  aov(data = meta,
      formula = n_perc_ww ~ treatment)
broom::tidy(mod_N) %>% 
  write_csv("./output/N_vs_treatment_aov-model.csv")
# P 
mod_P <- 
  aov(data = meta,
      formula = p_ppm ~ treatment)
broom::tidy(mod_P) %>% 
  write_csv("./output/P_vs_treatment_aov-model.csv")


## Model plots ####
# pH
ph_plot <- 
treatment_meta_plot("pH") +
  geom_label(data = multcomp_letters(mod_pH),aes(x=treatment,y=6,label=group),face='bold',size=5) +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# C
C_plot <- 
treatment_meta_plot("c_perc_ww") +
  geom_label(data = multcomp_letters(mod_C),aes(x=treatment,y=6,label=group),face='bold',size=5) +
  labs(y="C (% ww)") +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# H
H_plot <- 
treatment_meta_plot("h_perc_ww") +
  geom_label(data = multcomp_letters(mod_H),aes(x=treatment,y=2.0,label=group),face='bold',size=5) +
  labs(y="H (% ww)") +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# N
N_plot <- 
  treatment_meta_plot("n_perc_ww") +
  geom_label(data = multcomp_letters(mod_N),aes(x=treatment,y=0.9,label=group),face='bold',size=5) +
  labs(y="N (% ww)") +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# P
P_plot <- 
  treatment_meta_plot("p_ppm") +
  geom_label(data = multcomp_letters(mod_P),aes(x=treatment,y=800,label=group),face='bold',size=5) +
  labs(y="P (ppm)") +
  scale_x_discrete(labels=c('oil palm' = 'oil palm',
                            '0_1 years'='0-1 yr',
                            '1_2 years'='1-2 yr',
                            '2_3 years'='2-3 yr',
                            '3_4 years'='3-4 yr'))
ph_plot / C_plot / H_plot / N_plot / P_plot
ggsave("./output/figs/soil_vars_vs_treatment_w_stat-groups.png",dpi=400,height = 15,width = 8)  

# LEFTOVER CODE FOR STAT COMPARISONS IN A GIVEN PLOT
# p <- ggboxplot(meta, x = "treatment", y = "pH",
#                fill = "treatment", palette = pal$pal.okabe[c(2,4,6,7,8)],
#                add = "jitter") +
#   labs(x="\nTreatment") +
#   theme(legend.position = 'none',
#         line = element_line(linewidth = 6),
#         axis.title = element_text(face='bold',size=18),
#         axis.text = element_text(face='bold',size=14))
# p
# my_comparisons <- list( c("oil palm", "0-1 years"), 
#                         c("0-1 years", "1-2 years"), 
#                         c("2-3 years", "3-4 years"),
#                         c("oil palm", "1-2 years")
#                         )
# p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   stat_compare_means(label.y = 5,label.x = 3.6) +
#   scale_fill_viridis_d(option = 'Turbo',end=.8)
# ggsave("./output/figs/pH_vs_treatment.png",dpi=400,height = 6,width = 6)




# ALPHA-DIV VS METADATA ####

## by treatment ####
 
### fungi ####
# merge data and transform
fun_trt_ra <- 
fung %>%  
  merge_samples('treatment',fun = 'sum') %>% 
  transform_sample_counts(ra)
# repair relevant metadata
fix.fun <- 
  function(x){
    if(is.numeric(x)){
      return(x/phyloseq::nsamples(fun_trt_ra))} else {return(x)}
  }
fun_trt_ra@sam_data$treatment <- sample_names(fun_trt_ra)
fun_trt_ra@sam_data <- 
  map(fun_trt_ra@sam_data,fix.fun) %>% 
  as.data.frame() %>% 
  sample_data()
fun_trt_ra@sam_data$treatment <- 
  factor(fun_trt_ra@sam_data$treatment,levels = levels(meta$treatment) %>% str_replace('_','-'))

plot_bar2(fun_trt_ra,fill="Phylum",x='treatment') +
  scale_fill_viridis_d(option = 'turbo')
ggsave("./output/figs/fungi_phylum_by_treatment_barplot.png",dpi=400,height = 6,width = 5)

### bacteria ####

# merge data and transform
bac_trt_ra <- 
  bact %>%  
  merge_samples('treatment',fun = 'sum') %>% 
  transform_sample_counts(ra)
# repair relevant metadata
fix.fun <- 
  function(x){
    if(is.numeric(x)){
      return(x/phyloseq::nsamples(bac_trt_ra))} else {return(x)}
  }
bac_trt_ra@sam_data$treatment <- sample_names(bac_trt_ra)
bac_trt_ra@sam_data <- 
  map(bac_trt_ra@sam_data,fix.fun) %>% 
  as.data.frame() %>% 
  sample_data()
bac_trt_ra@sam_data$treatment <- 
  factor(bac_trt_ra@sam_data$treatment,levels = levels(meta$treatment) %>% str_replace('_','-'))

plot_bar2(bac_trt_ra,fill="Phylum",x='treatment') +
  scale_fill_viridis_d(option = 'turbo')
ggsave("./output/figs/bacteria_phylum_by_treatment_barplot.png",dpi=400,height = 6, width = 10)

# plot richness by group
richness <- 
full_join(
  bact %>% 
    estimate_richness(measures = "Observed") %>% 
    mutate(treatment = fung@sam_data$treatment,
           taxa="Bacteria"),
  fung %>% 
    estimate_richness(measures = "Observed") %>% 
    mutate(treatment = fung@sam_data$treatment,
           taxa="Fungi")
)
richness %>% 
  ggplot(aes(x=treatment,y=Observed,fill=treatment)) +
  geom_boxplot(color='black') +
  geom_jitter(width = .05) + 
  facet_wrap(~taxa,scales='free') +
  scale_fill_manual(values = trt.pal) +
  labs(fill="Treatment",x="\nTreatment") +
  theme(strip.background = element_rect(fill='white'),
        axis.text.x = element_text(angle = 90,hjust=1,vjust = .5))
ggsave("./output/figs/Richness_by_treatment_boxplots.png",dpi=400,height = 8,width = 12)




## by soil parameters ####
# (use plot as blocking, random effect)
# rescale predictors
scaled_meta <- 
  meta %>%
  mutate(across(all_of(c("pH" , "c_perc_ww" , "h_perc_ww" , "n_perc_ww" , "p_ppm")),scale))

# fungal ASV richness
fung_asv_rich_mod <- 
  lmerTest::lmer(data=scaled_meta,
                 formula=fungal_asv_richness ~ pH + c_perc_ww + h_perc_ww + n_perc_ww + p_ppm + (1|plot))
lmerTest:::get_coefmat(fung_asv_rich_mod) %>% 
  as.data.frame() %>% 
  write_csv("./output/fungal_asv_richness_vs_soil_lmer-mod.csv")
# bacterial ASV richness
bact_asv_rich_mod <- 
  lmerTest::lmer(data=scaled_meta,
                 formula=bacterial_asv_richness ~ pH + c_perc_ww + h_perc_ww + n_perc_ww + p_ppm + (1|plot))
lmerTest:::get_coefmat(bact_asv_rich_mod) %>% 
  as.data.frame() %>% 
  write_csv("./output/bacterial_asv_richness_vs_soil_lmer-mod.csv")

# fungal spp richness
fung_spp_rich_mod <- 
  lmerTest::lmer(data=scaled_meta,
                 formula=fungal_spp_richness ~ pH + c_perc_ww + h_perc_ww + n_perc_ww + p_ppm + (1|plot))
lmerTest:::get_coefmat(fung_spp_rich_mod) %>% 
  as.data.frame() %>% 
  write_csv("./output/fungal_spp_richness_vs_soil_lmer-mod.csv")
# bacterial spp richness
bact_spp_rich_mod <- 
  lmerTest::lmer(data=scaled_meta,
                 formula=bacterial_spp_richness ~ pH + c_perc_ww + h_perc_ww + n_perc_ww + p_ppm + (1|plot))
lmerTest:::get_coefmat(bact_spp_rich_mod) %>% 
  as.data.frame() %>% 
  write_csv("./output/bacterial_spp_richness_vs_soil_lmer-mod.csv")

### plots ####
meta %>% 
  pivot_longer(ends_with('asv_richness'),
               names_to = 'Domain',values_to = "Richness") %>% 
  mutate(Domain = case_when(Domain == 'fungal_asv_richness' ~ "Fungi",
                            Domain == 'bacterial_asv_richness' ~ "Bacteria")) %>% 
  pivot_longer(all_of(c('pH','c_perc_ww','h_perc_ww','n_perc_ww','p_ppm')),
               names_to = "Variable",values_to = "Value") %>% 
  mutate(Variable = Variable %>% 
           str_to_sentence() %>% 
           str_replace("_perc_ww"," (% ww)") %>% 
           str_replace("_ppm"," (ppm)") %>% 
           str_replace("Ph","pH")) %>% 
  ggplot(aes(x=Value,y=Richness,color=Domain)) +
  geom_point(size=3,alpha=.5) +
  geom_smooth(se=FALSE,method='lm',linewidth=2) +
  facet_grid(vars(Domain),vars(Variable),scales = 'free') +
  scale_color_manual(values = pal$pal.earthtones) +
  labs(y="ASV richness")
ggsave("./output/figs/ASV_Diversity_both-domains_vs_soil_vars.png",dpi=400,height = 6,width = 12)

meta %>% 
  pivot_longer(ends_with('spp_richness'),
               names_to = 'Domain',values_to = "Richness") %>% 
  mutate(Domain = case_when(Domain == 'fungal_spp_richness' ~ "Fungi",
                            Domain == 'bacterial_spp_richness' ~ "Bacteria")) %>% 
  pivot_longer(all_of(c('pH','c_perc_ww','h_perc_ww','n_perc_ww','p_ppm')),
               names_to = "Variable",values_to = "Value") %>% 
  mutate(Variable = Variable %>% 
           str_to_sentence() %>% 
           str_replace("_perc_ww"," (% ww)") %>% 
           str_replace("_ppm"," (ppm)") %>% 
           str_replace("Ph","pH")) %>% 
  ggplot(aes(x=Value,y=Richness,color=Domain)) +
  geom_point(size=3,alpha=.5) +
  geom_smooth(se=FALSE,method='lm',linewidth=2) +
  facet_grid(vars(Domain),vars(Variable),scales = 'free') +
  scale_color_manual(values = pal$pal.earthtones) +
  labs(y="Species richness")
ggsave("./output/figs/Spp_Diversity_both-domains_vs_soil_vars.png",dpi=400,height = 6,width = 12)

full_join(
  meta %>% 
    pivot_longer(ends_with('spp_richness'),
                 names_to = 'Domain',values_to = "Richness") %>% 
    mutate(Domain = case_when(Domain == 'fungal_spp_richness' ~ "Fungi",
                              Domain == 'bacterial_spp_richness' ~ "Bacteria")) %>% 
    pivot_longer(all_of(c('pH','c_perc_ww','h_perc_ww','n_perc_ww','p_ppm')),
                 names_to = "Variable",values_to = "Value") %>% 
    mutate(Variable = Variable %>% 
             str_to_sentence() %>% 
             str_replace("_perc_ww"," (% ww)") %>% 
             str_replace("_ppm"," (ppm)") %>% 
             str_replace("Ph","pH")) %>%
    dplyr::filter(Domain == "Fungi") %>% 
    glm(data=.,formula=Richness ~ Value * Variable) %>% 
    broom::tidy() %>% 
    mutate(Domain = "Fungi"),
  meta %>% 
    pivot_longer(ends_with('spp_richness'),
                 names_to = 'Domain',values_to = "Richness") %>% 
    mutate(Domain = case_when(Domain == 'fungal_spp_richness' ~ "Fungi",
                              Domain == 'bacterial_spp_richness' ~ "Bacteria")) %>% 
    pivot_longer(all_of(c('pH','c_perc_ww','h_perc_ww','n_perc_ww','p_ppm')),
                 names_to = "Variable",values_to = "Value") %>% 
    mutate(Variable = Variable %>% 
             str_to_sentence() %>% 
             str_replace("_perc_ww"," (% ww)") %>% 
             str_replace("_ppm"," (ppm)") %>% 
             str_replace("Ph","pH")) %>%
    dplyr::filter(Domain == "Bacteria") %>% 
    glm(data=.,formula=Richness ~ Value * Variable) %>% 
    broom::tidy() %>% 
    mutate(Domain = "Bacteria")
) %>% 
  mutate(term = term %>% str_remove("Variable")) %>% 
  write_csv("output/Spp_richness_glm_mod_table_full.csv")


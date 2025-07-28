# BETA DIVERSITY 

# SETUP ####
set.seed(666)

## packages ####
library(tidyverse)
library(phyloseq)
library(janitor)
library(patchwork)
library(broom)
library(easystats)
library(vegan)
library(rbiom)
library(ape)
library(TreeTools)
library(gdm)


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


## load physeq objects ####
fung <- readRDS("./output/ITS_Physeq_cleaned_w_tree.RDS")
bact <- readRDS("./output/16S_Physeq_cleaned_w_tree.RDS")

# ORDINATIONS ####
fung_ord <- fung %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  ordinate(method="NMDS",distance = "unifrac")
bact_ord <- bact %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  ordinate(method="NMDS",distance = "unifrac")
fung_ord$stress
bact_ord$stress

bact_nmds <- plot_ordination(bact, bact_ord, color="treatment") + 
  geom_point(size=6,alpha=.8) +
  scale_color_manual(values=trt.pal) +
  scale_shape_manual(values = c(15, 16, 17, 18, 19)) +
  labs(color="Treatment",title = "Bacteria") +
  theme(legend.position = 'none') +
  stat_ellipse(level = .75,linetype=221)

fung_nmds <- plot_ordination(fung, fung_ord, color="treatment",shape='treatment') + 
  geom_point(size=6,alpha=.8) +
  scale_color_manual(values=trt.pal) +
  scale_shape_manual(values = c(15, 16, 17, 18, 19)) +
  labs(color="Treatment",shape="Treatment",title = "Fungi") +
  stat_ellipse(level = .75,linetype=221)

# joint plot  
bact_nmds | fung_nmds
ggsave("./output/figs/NMDS_Plot_by_treatment.png",dpi=400,height = 8,width = 16)

# PERMANOVA ####
adonis2(otu_table(fung) ~ fung@sam_data$treatment) %>% 
  tidy() %>% 
  mutate(term=term %>% str_remove('fung@sam_data\\$'),domain="Fungi") %>% 
  full_join(
    adonis2(otu_table(bact) ~ bact@sam_data$treatment) %>% 
      tidy() %>% 
      mutate(term=term %>% str_remove('bact@sam_data\\$'),domain="Bacteria")
  ) %>% 
  write_csv("./output/permanova_table_treatment.csv")

adonis2(otu_table(tax_glom(fung, taxrank = "Species")) ~ fung@sam_data$treatment) %>% 
  tidy() %>% 
  mutate(term=term %>% str_remove('fung@sam_data\\$'),domain="Fungi") %>% 
  full_join(
    adonis2(otu_table(tax_glom(bact, taxrank = "Species")) ~ bact@sam_data$treatment) %>% 
      tidy() %>% 
      mutate(term=term %>% str_remove('bact@sam_data\\$'),domain="Bacteria")
  ) %>% 
  write_csv("./output/permanova_table_treatment_spp-level.csv")

# BETA-DISPERSION ####

dist_fung <- rbiom::unifrac(biom=fung %>% transform_sample_counts(ra) %>% otu_table() %>% as('matrix') %>% t(),
                            weighted = TRUE,
                            tree=fung@phy_tree)
dist_bact <- rbiom::unifrac(biom=bact %>% transform_sample_counts(ra) %>% otu_table() %>% as('matrix') %>% t(),
                            weighted = TRUE,
                            tree=bact@phy_tree)
bd_fung <- betadisper(d = dist_fung,
           group = fung@sam_data$treatment)
bd_bact <- betadisper(d = dist_bact,
                      group = bact@sam_data$treatment)
plot(bd_bact);plot(bd_fung)

# make complete data frame of dispersion values
beta_disp_df <- 
data.frame(
  Group = c(bd_fung$group,bd_bact$group),
  Dispersion = c(bd_fung$distances,bd_bact$distances),
  Domain = c(rep("Fungi",length(bd_fung$distances)),
             rep("Bacteria",length(bd_bact$distances)))
)
# plot beta-disp
beta_disp_df %>% 
  ggplot(aes(x=Group,y=Dispersion,fill=Domain)) +
  geom_boxplot() +
  scale_fill_manual(values = pal$pal.earthtones) +
  labs(x="\nTreatment",y="Beta-dispersion") +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5))
ggsave("./output/figs/Beta-dispersion_plot.png",height = 5,width = 8,dpi=400)
# model beta-disp by location
bact_beta.disp_mod <- 
  glm(data=beta_disp_df %>% filter(Domain == "Bacteria"),
      formula=Dispersion ~ Group)
fung_beta.disp_mod <- 
  glm(data=beta_disp_df %>% filter(Domain == "Fungi"),
      formula=Dispersion ~ Group)


broom::tidy(bact_beta.disp_mod) %>% 
  mutate(Domain="Bacteria") %>% 
  full_join(
    broom::tidy(fung_beta.disp_mod) %>% 
      mutate(Domain="Fungi")
  ) %>% mutate(term=term %>% str_remove("Group")) %>% 
  write_csv("./output/Beta-dispersion_model.csv")




# GDM MODEL ####
# extract species by site info
fung_ra_melt <- 
  fung %>% 
  transform_sample_counts(ra) %>% 
  psmelt()
bact_ra_melt <- 
  bact %>% 
  transform_sample_counts(ra) %>% 
  psmelt()


# biological data
# get columns with xy, site ID, and species data
sppTab_fung <- fung_ra_melt %>% dplyr::select(OTU,Sample,longitude,latitude,treatment,Abundance)
sppTab_bact <- bact_ra_melt %>% dplyr::select(OTU,Sample,longitude,latitude,treatment,Abundance)
# get columns with site ID, env. data, and xy-coordinates
envTab_fung <- fung_ra_melt %>% dplyr::select(Sample,longitude,latitude,pH,ends_with("perc_ww"),ends_with("_ppm"))
envTab_bact <- bact_ra_melt %>% dplyr::select(Sample,longitude,latitude,pH,ends_with("perc_ww"),ends_with("_ppm"))
# format for gdm
gdmTab_fung <- formatsitepair(bioData=sppTab_fung, weightType = 'richness',
                             bioFormat=2, #x-y spp list
                             XColumn="longitude", 
                             YColumn="latitude",
                             sppColumn="OTU", 
                             siteColumn="Sample", 
                             predData=envTab_fung,
                             abundance = TRUE,
                             abundColumn = "Abundance")
gdmTab_bact <- formatsitepair(bioData=sppTab_bact, weightType = 'richness',
                              bioFormat=2, #x-y spp list
                              XColumn="longitude", 
                              YColumn="latitude",
                              sppColumn="OTU", 
                              siteColumn="Sample", 
                              predData=envTab_bact,
                              abundance = TRUE,
                              abundColumn = "Abundance")

# fit GDM models
gdm_fung <- gdm(data = gdmTab_fung,geo = TRUE)
gdm_bact <- gdm(data = gdmTab_bact,geo = TRUE)
# quick look at model fits
summary(gdm_fung)
summary(gdm_bact)

# predictions from model (using same distances)
gdm_fung_pred <- predict(object=gdm_fung, data=gdmTab_fung)
gdm_bact_pred <- predict(object=gdm_bact, data=gdmTab_bact)

fung_preds <- data.frame(observed = gdmTab_fung$distance,
                         predicted = gdm_fung_pred,
                         dist = gdm_fung$ecological,
                         sample_type = "Fungi")
bact_preds <- data.frame(observed = gdmTab_bact$distance,
                         predicted = gdm_bact_pred,
                         dist = gdm_bact$ecological,
                         sample_type = "Bacteria")
full_join(fung_preds,bact_preds) %>% 
  ggplot(aes(x=observed,y=predicted)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_wrap(~sample_type,scales = 'free')
ggsave("./output/figs/gdm_obs_vs_predicted.png",width = 8,height = 6,dpi=300)
full_join(fung_preds,bact_preds) %>% 
  ggplot(aes(x=dist,y=observed)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~sample_type,scales = 'free')
ggsave("./output/figs/gdm_comm-dist_vs_observed.png",width = 8,height = 6,dpi=300)

## Extract splines ####
fung_splines <- gdm::isplineExtract(gdm_fung) %>% as.data.frame() %>% mutate(Domain='Fungi')
bact_splines <- gdm::isplineExtract(gdm_bact) %>% as.data.frame() %>% mutate(Domain='Bacteria')

isplines <- full_join(bact_splines,fung_splines)



names(isplines)
names(isplines) <- c("geographic_actual","pH_actual","C_actual","H_actual","N_actual","P_actual",
                     "geographic_partial","pH_partial","C_partial","H_partial","N_partial","P_partial","Domain")
## Plot splines ####
p1 <- isplines %>% 
  ggplot(aes(x=pH_actual,y=pH_partial,color=Domain)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2) +
  scale_color_manual(values=pal$pal.earthtones) +
  labs(y="f(pH)",x="pH",color="Taxa")
p1

p2 <- isplines %>% 
  ggplot(aes(x=geographic_actual*10000,y=geographic_partial,color=Domain)) +
  geom_smooth(se=FALSE,linewidth=2) +
  scale_color_manual(values=pal$pal.earthtones) +
  labs(y="f(Geo. dist.)",x="Geographic distance (m)",color="Taxa")

p3 <- isplines %>% 
  ggplot(aes(x=C_actual,y=C_partial,color=Domain)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2) +
  scale_color_manual(values=pal$pal.earthtones) +
  labs(y="f(C)",x="C (% ww)",color="Taxa")

p4 <- isplines %>% 
  ggplot(aes(x=H_actual,y=H_partial,color=Domain)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2) +
  scale_color_manual(values=pal$pal.earthtones) +
  labs(y="f(H)",x="H (% ww)",color="Taxa")

p5 <- isplines %>% 
  ggplot(aes(x=N_actual,y=N_partial,color=Domain)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2) +
  scale_color_manual(values=pal$pal.earthtones) +
  labs(y="f(N)",x="N (% ww)",color="Taxa")

p6 <- isplines %>% 
  ggplot(aes(x=P_actual,y=P_partial,color=Domain)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2) +
  scale_color_manual(values=pal$pal.earthtones) +
  labs(y="f(P)",x="P (ppm)",color="Taxa")

(p2 | p1) / (p3 | p4) / (p5 | p6) + plot_layout(guides='collect')
ggsave("./output/figs/GDM_splines_plots_soilvars.png",dpi=400,height = 8,width = 12)




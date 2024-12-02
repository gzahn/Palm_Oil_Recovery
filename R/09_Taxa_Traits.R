# COMMUNITY ASSEMBLY OVER TIME

# SETUP ####

## Packages ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(patchwork); packageVersion("patchwork")
library(BacDive); packageVersion("BacDive")
library(fungaltraits); packageVersion("fungaltraits")

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

# LOAD DATA ####
# Load cleaned phyloseq object
ps <- readRDS("./output/clean_phyloseq_object.RDS")

# load list of significant taxa
# sig_taxa <- readRDS("./output/final_significant_taxa.RDS")

                    ############################################################
                    # Stuff below was for a different project...needs adapting #
                    ############################################################

# Fungal Traits ####

# download traits metadata
traits_meta <- read_csv("https://github.com/traitecoevo/fungaltraits/releases/download/v0.0.3/funtothefun.csv")

# download FungalTraits database
traits_db <- fungaltraits::fungal_traits()
names(traits_db$species)
# match taxa at genus level
genera <- fung@tax_table[,6] %>% str_remove("^g__")
species <- fung@tax_table[,7] %>% str_remove("^s__")
fungal_traits <- 
  data.frame(Genus=genera) %>% 
  mutate(species=paste(Genus,species,sep="_")) %>% 
  left_join(traits_db,by=c("species","Genus"),multiple='all')

# need to condense/remove multiple matches
fungal_traits %>% 
  dplyr::filter(species != "NA_NA")

# remove traits not associated with biochem functional potential
traits_to_ignore <- c(
  "redChannel_mean","redChannel_sd","RNAHelicase_count","RNApolymerase_count","spore_length",
  "spore_size","spore_width","sporocarp_chitin","sporocarp_N","sporocarp_protein","sporocarp_resp",           
  "taxonomic_level_fg","tissue_c","tissue_cn","tissue_cp","tissue_n","tissue_np","tissue_p","total_genes",
  "trehalase_count","latitude","map","greenChannel_mean","greenChannel_sd","heatShockProtein_count",
  "extension_rate","fruiting_body_size","mat","longitude","melanin_content","melanin_count",
  "coldShockProtein_count","dsDNA","blueChannel_mean","blueChannel_sd","ifungorum_number",
  "sterol_type","studyName","substrate","trait_fg","trophic_mode_fg",'notes_fg',"source_funguild_fg",
  "growth_form_fg","guild_fg","higher_clade","culture_media","culture_notes","elevation","em_expl",
  "em_text","colour_mean","confidence_fg","ascoma_development","ascoma_type","ascus_dehiscence",
  "uuid","obj_id","speciesMatched"
)

# group by species; summarize to find mean values with na.omit=TRUE
summarized_traits <- 
  fungal_traits %>% 
  dplyr::select(-all_of(traits_to_ignore)) %>% 
  dplyr::group_by(species) %>% 
  summarize(across(where(is.numeric),function(x){mean(x,na.rm=TRUE)}))

names(summarized_traits)

# join traits with tax_table species 

traits <- 
  data.frame(Genus=genera) %>% 
  mutate(species=paste(Genus,species,sep="_")) %>% 
  left_join(summarized_traits,by=c("species"))




# BACDIVE ####
# initialize bacdive connection (add your userid and password)
# BacDive URL: https://bacdive.dsmz.de/
bacdive <- BacDive::open_bacdive(username = Sys.getenv("BACDIVE_USER"),
                                 password = Sys.getenv("BACDIVE_PW"))



# get list of unique genera
genus_list <- ps@tax_table[,6] %>% 
  table()
genus_list <- genus_list %>% 
  as.data.frame() %>% 
  pluck(".") %>% 
  levels() 

# BacDive API ####

# build large list of BacDive cell morphology
morph_list <- list()

for(i in genus_list){
  df <- get_bacdive_morphology(i)
  morph_list[[i]] <- df
}  
saveRDS(morph_list,"./output/genus_morphology_data.RDS")
morph_list <- readRDS("./output/genus_morphology_data.RDS")

names(morph_list)
# morph_list <- readRDS("./output/genus_morphology_data.RDS")

# reduce to single data frame
morphology <- purrr::reduce(morph_list,full_join)

# ANALYSIS OF BACDIVE DATA ####
# minimum dimensions
morphology$min_length <- morphology$length %>% 
  str_remove(" µm") %>% 
  str_split("-") %>% 
  map_chr(1) %>% 
  as.numeric()

morphology$min_width <- morphology$width %>% 
  str_remove(" µm") %>% 
  str_split("-") %>% 
  map_chr(1) %>% 
  as.numeric()


# add signifigance column
morphology <- morphology %>% 
  mutate(signifigant_taxa = case_when(genus %in% sig_taxa ~ TRUE,
                                      TRUE ~ FALSE),
         min_dimension = ifelse((min_length - min_width) > 0, 
                                min_width,
                                min_length),
         est_vol = min_length * min_width,
         shape = shape %>% str_remove("-shaped"),
         est_surface_area = case_when(shape == "coccus" ~ 4*pi*((min_dimension/2)^2),
                                      shape != "coccus" ~ (2*pi*((min_length/2)^2))))

saveRDS(morphology,"./output/cell_morphology_data.RDS")
ps@sam_data
morphology %>% 
  dplyr::filter(signifigant_taxa) %>% 
  group_by(genus) %>% 
  summarize(Avg_Min = mean(min_length,na.rm=TRUE))
write_csv("./output/significant_taxa_traits.csv")


morphology %>% 
  dplyr::filter(signifigant_taxa) %>% 
  mutate(length=length %>% str_remove(" µm") %>% str_split("-") %>% map_chr(1) %>% as.numeric,
         width=width %>% str_remove(" µm") %>% str_split("-") %>% map_chr(1) %>% as.numeric) %>% 
  group_by(genus) %>% 
  summarize(mean_length=mean(length,na.rm=TRUE),
            mean_width=mean(width,na.rm=TRUE))

# distribution plots
morphology %>% 
  ggplot(aes(x=min_dimension, fill = signifigant_taxa)) +
  geom_density(alpha=.5) +
  labs(x="Minimum dimension (µm)",
       fill="Significant\ntaxa")
ggsave("./output/figs/minimum_cell_dimension_distribution.png",
       dpi=300,height = 6, width = 6)

# GLM model output for minimum cell dimension
# minimum dimension is important for Reynold's Number
mod <- glm(data=morphology,
           formula = min_dimension ~ signifigant_taxa)
saveRDS(mod,"./output/cell_dimension_glm.RDS")
sink("./output/cell_min_dimension_glm_summary.txt")
mod %>% summary()
sink(NULL)
report::report(mod)
morphology
# Chi-Square test for enrichment in cell shape for significant taxa

xsqtest <- table(morphology$shape, morphology$signifigant_taxa) %>% 
  chisq.test()
saveRDS(xsqtest,"./output/cell_shape_xsq.RDS")

sink("./output/cell_shape_chi-sq_test.txt")
print("Genus significance and shape")
table(morphology$shape, morphology$signifigant_taxa)
print("")
table(morphology$shape, morphology$signifigant_taxa) %>% 
  chisq.test()
print("No apparent enrichment in cell shape morphology in significant taxa")
sink(NULL)

morphology$est_surface_area
table(morphology$shape)


# surface area 
morphology %>% 
  ggplot(aes(x=est_surface_area,fill=signifigant_taxa)) +
  geom_density()

mod_surfacearea <- 
  glm(data=morphology,
      formula = est_surface_area ~ signifigant_taxa)
report::report(mod_surfacearea)
mod_surfacearea %>% summary
sink("./output/cell_surface_area_glm_summary.txt")
summary(mod_surfacearea)
sink(NULL)


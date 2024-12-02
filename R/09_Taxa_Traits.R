# COMMUNITY ASSEMBLY OVER TIME

# SETUP ####

## Packages ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(patchwork); packageVersion("patchwork")
library(BacDive); packageVersion("BacDive")
library(fungaltraits); packageVersion("fungaltraits")
library(FUNGuildR); packageVersion("FUNGuildR")
library(lemon); packageVersion("lemon")

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
bact <- readRDS("./output/16S_Physeq_cleaned_w_tree.RDS")
fung <- readRDS("./output/ITS_Physeq_cleaned_w_tree.RDS")

# load list of significant taxa
hub.taxa <- read_csv("./output/hub_taxa_list.csv")
new.arrivals <- read_csv("./output/new_arrival_taxa_over_time.csv")
fung.sig.taxa <- readRDS("./output/fungi_significant_taxa.RDS")
bact.sig.taxa <- readRDS("./output/bacteria_significant_taxa.RDS")                    


                    ############################################################
                    # Stuff below was for a different project...needs adapting #
                    ############################################################

# Fungal Traits ####

## FunGuild ##
guild_db <- FUNGuildR::get_funguild_db()

# save guild db as RDS
saveRDS(guild_db, "./taxonomy/Funguild_Database.RDS")

# assign guild to fungal ASV taxonomy
guilds <- 
  data.frame(
    Taxonomy = paste(
      tax_table(fung)[,1],
      tax_table(fung)[,2],
      tax_table(fung)[,3],
      tax_table(fung)[,4],
      tax_table(fung)[,5],
      tax_table(fung)[,6],
      tax_table(fung)[,7],
      sep=";"
    )
  ) %>% 
  FUNGuildR::funguild_assign()

guilds
# add guild info and rename ASVs to match
fung@tax_table[,1] <- guilds$guild
taxa_names(fung) <- paste0("ASV_",seq_along(taxa_names(fung)))
taxa_names(bact) <- paste0("ASV_",seq_along(taxa_names(bact)))


fung.sig.melt <- 
  fung %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(taxa_names(fung) %in% fung.sig.taxa$OTU) %>% 
  psmelt()
fung.sig.melt <- 
  fung.sig.melt %>%
  mutate(Guild = case_when(grepl("Plant Pathogen",Kingdom) ~ "Plant Pathogen",
                           grepl("Fungal Parasite",Kingdom) ~ "Fungal Parasite",
                           grepl("Animal Pathogen",Kingdom) ~ "Animal Pathogen",
                           grepl("Saprotroph",Kingdom) ~ "Saprotroph"))
  
# remake important fungi plots
p <- 
fung.sig.melt %>% 
  dplyr::filter(!is.na(Genus) & Genus != "unclassified" & !is.na(Guild)) %>% 
  # dplyr::filter(treatment != "oil palm") %>% 
  ggplot(aes(x=treatment,y=Abundance,fill=Guild)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5),
        strip.background = element_rect(fill='white'),
        strip.text = element_text(size=12,face='bold.italic'),
        panel.spacing = unit(1,"lines"),
        strip.placement = 'outside'
       ) +
  scale_fill_viridis_d(option="turbo") +
  labs(x="Time since restoration",y="Relative abundance") +
  facet_wrap(~Genus,scales = 'free_y',axis.labels = 'margins')
  
p  

ggsave("./output/figs/important_fungal_taxa_guilds.png",height = 16, width = 16, dpi=400)
# these fungi increased since oil palm, at least in one year



# BACDIVE ####
# initialize bacdive connection (add your userid and password)
# BacDive URL: https://bacdive.dsmz.de/
bacdive <- BacDive::open_bacdive(username = Sys.getenv("BACDIVE_USER"),
                                 password = Sys.getenv("BACDIVE_PW"))

bact.sig.melt <- 
  bact %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(taxa_names(bact) %in% bact.sig.taxa$OTU) %>% 
  psmelt()


p <- 
  bact.sig.melt %>% 
  dplyr::filter(!is.na(Genus) & Genus != "unclassified" ) %>% 
  # dplyr::filter(treatment != "oil palm") %>% 
  ggplot(aes(x=treatment,y=Abundance,fill=Phylum)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5),
        strip.background = element_rect(fill='white'),
        strip.text = element_text(size=12,face='bold.italic'),
        panel.spacing = unit(1,"lines"),
        strip.placement = 'outside'
  ) +
  scale_fill_viridis_d(option="turbo") +
  labs(x="Time since restoration",y="Relative abundance") +
  facet_wrap(~Genus,scales = 'free_y',axis.labels = 'margins')

p  

ggsave("./output/figs/important_bacterial_taxa_guilds.png",height = 16, width = 16, dpi=400)
# get list of unique genera to look up in literature
genus_list <- bact.sig.melt$Genus %>% str_remove("Candidatus ") %>% unique
genus_list <- genus_list[!is.na(genus_list)]
writeLines(genus_list,"./output/significant_bacterial_genera.txt")
genus_list

# Udaeobacter is one of the most cosmopolitan soil fungi. Genetic streamlining.
# Halangium produce the antifungal compounds haliangicins.
# Occallatibacter acidobacteria isolated from Namibian soils
# Bacillus strains here are bacillus living in salty or alkaline conditions, or salty and alkaline conditions
#   "megaterium" can solubilize phosphate
# ADurb.Bin063-1 has been noted for Annamox ability
# Koribacter CO-oxidizing, potential large role in CO removal from atmosphere
# Reyranella ???
# Aquisphaera optimum growth temperature of about 30-35 °C and an optimum pH for growth of around 7.5-8.5
# 

x <- readRDS("./output/16S_Physeq_cleaned_w_tree.RDS")
# x[,7] %>% unique %>% unname
plot_P_taxa <- 
function(m,g,s){
  t <- try(m %>% subset_taxa(Genus==g & Species==s) %>% ntaxa)
  if("try-error" %in% class(t)){stop("Taxon not in data set.")} 
  
  if("try-error" %ni% class(t)) {
    melt <- psmelt(m %>% transform_sample_counts(ra))
    x <- 
      melt %>% 
      dplyr::filter(!is.na(Genus) & Genus != "unclassified") %>%
      dplyr::filter(Genus == g & Species == s) 
    
    x$SciName <- paste(x$Genus,x$Species)  
    x %>% 
      ggplot(aes(x=treatment,y=Abundance)) +
      geom_boxplot() +
      theme(axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5),
            strip.background = element_rect(fill='white'),
            strip.text = element_text(size=12,face='bold.italic'),
            panel.spacing = unit(1,"lines"),
            strip.placement = 'outside'
      ) +
      scale_fill_viridis_d(option="turbo") +
      labs(x="Time since restoration",y="Relative abundance") +
      facet_wrap(~SciName,scales = 'free_y',axis.labels = 'margins')
    ggsave(paste0("./output/figs/Sig_Taxa_P_solubilizers_",g,"_",s,".png"),height = 6,width = 6,dpi=300)
    
  }
  
}

# m=bact;g="Bacillus";s="circulans"
plot_P_taxa(bact,"Bacillus","circulans")
plot_P_taxa(bact,"Bacillus","megaterium")
plot_P_taxa(bact,"Bacillus","megaterium")

# do with fungi as well
f <- readRDS("./output/ITS_Physeq_cleaned_w_tree.RDS")
plot_P_taxa(f,"Aspergillus","fumigatus")
plot_P_taxa(f,"Trichoderma","viride")


bm <- bact %>% 
  transform_sample_counts(ra) %>%
  subset_taxa(Genus == "Bacillus" & Species == "megaterium") %>% 
  psmelt()

bm %>% 
  ggplot(aes(x=Abundance,y=p_ppm)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x="Relative abundance",y="P (ppm)",title = "Bacillus megaterium") +
  theme(plot.title = element_text(face='bold.italic',size=14,hjust=.5))
ggsave("./output/figs/Bacillus_megaterium_P_solubilizer_vs_Pppm.png",height = 6,width = 6,dpi=300)
glm(data=bm,
    formula = p_ppm ~ Abundance) %>% broom::tidy() %>% 
  write_csv("./output/figs/Bacillus_megaterium_P_solubilizer_vs_Pppm_glm.csv")

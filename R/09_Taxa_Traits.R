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


    


# FUNGUILD ####
guild_db <- FUNGuildR::get_funguild_db()

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
    ),
    ASV = taxa_names(fung)
  ) %>% 
  FUNGuildR::funguild_assign() %>% 
  mutate(trophicMode = ifelse(is.na(trophicMode),"Undefined",trophicMode))

guilds %>% 
  saveRDS("./output/fungal_trait_df.RDS")
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

# ggsave("./output/figs/important_fungal_taxa_guilds.png",height = 16, width = 16, dpi=400)
# these fungi increased since oil palm, at least in one year



# BACDIVE ####
# initialize bacdive connection (add your userid and password)
# BacDive URL: https://bacdive.dsmz.de/
bacdive <- BacDive::open_bacdive(username = Sys.getenv("BACDIVE_USER"),
                                 password = Sys.getenv("BACDIVE_PW"))

genus <- "Achromobacter"
# bacdive helper function
get_bacdive_morphology <- function(genus){
  x <- BacDive::request(object = bacdive,
                        query = genus,
                        search = "taxon")
  # find bacdive id
  y <- x$results
  
  if(length(y) == 0){
    df <- data.frame(genus=genus,
                     known_enzymes=NA,
                     num_known_enzymes=NA)
    return(df)
  }
  
  # get info on that id
  z <- BacDive::fetch(bacdive,y)
  
  # get known enzyme list
  m <- z$results %>% 
    map("Physiology and metabolism") %>% 
    map("enzymes")
  
  nulls <- 
    m %>% 
    map(is.null) %>% 
    unlist
  m <- m[!nulls]
  m.values <- 
    m %>% 
    unlist
  known_enzymes <- 
    m.values[grep(pattern = "value$",names(m.values),value = TRUE)] %>% 
    unique

  if(length(known_enzymes) < 1){
    num_known_enzymes <- NA
  } else {
    num_known_enzymes <- length(known_enzymes)
  }
  
    
    df <- data.frame(genus=genus,
                     known_enzymes=paste(known_enzymes,collapse = ";"),
                     num_known_enzymes=num_known_enzymes)
    return(df)
  }
  




## Extract genus names ####
# get list of unique genera
genus_list <- bact@tax_table[,6] %>% 
  table()
genus_list <- genus_list %>% 
  as.data.frame() %>% 
  pluck(".") %>% 
  levels() 


## Query BacDive ####
# build large list of BacDive cell morphology
enzyme_list <- list()

for(i in genus_list){
  df <- get_bacdive_morphology(i)
  enzyme_list[[i]] <- df
}  

saveRDS(enzyme_list,"./output/bact_genus_enzyme_data.RDS")
genus_enzyme_df <- 
data.frame(genus=names(enzyme_list),
           num_known_enzymes=enzyme_list %>% map_dbl("num_known_enzymes"))
saveRDS(genus_enzyme_df,"./output/bacterial_enzyme_diversity_df.RDS")


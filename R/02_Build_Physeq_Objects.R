# Build Phyloseq Objects

# SETUP ####
set.seed(666)
## packages ####

# for HPC
mylib <- "/uufs/chpc.utah.edu/common/home/u6033249/R/library-4.4"

library(tidyverse, lib.loc = mylib)
library(dada2, lib.loc = mylib)
library(phyloseq, lib.loc = mylib)
library(decontam, lib.loc =  mylib)
library(janitor, lib.loc = mylib)

## functions ####
source("./R/functions.R")

## load metadata ####
meta <- read_csv("./data/Metadata.csv") %>% 
  janitor::clean_names()

## data checks ####
if(any(!file.exists(meta$filepaths))){
  warning("Some files not present in given paths. Check download. Subsetting to only present files.")
  meta <- meta[file.exists(meta$filepaths),]
}  

if(!file.exists("./taxonomy/Eukaryome_General_ITS_v1.8_reformatted.fasta.gz") |
   !file.exists("./taxonomy/Eukaryome_General_Longread_v1.9_reformatted.fasta.gz")){
  stop("Taxonomic databases not downloaded into correct path.")
}





# BUILD ASV TABLES ####

# for 16S
asv_tab_16S <- 
build_asv_table(metadata = meta,
                filtered_path = "./data/16S/filtered",
                amplicon = "16S",
                multithread = parallel::detectCores()-1,
                control.column = "control",
                neg.ctl.term = "Neg",
                out.dir = "./output")
# reduce asv table to reasonable asvs
asv_tab_16S <- asv_tab_16S[,colSums(asv_tab_16S) > 99]

# for ITS
asv_tab_ITS <- 
build_asv_table(metadata = meta,
                filtered_path = "./data/16S/filtered",
                amplicon = "ITS",
                multithread = parallel::detectCores()-1,
                control.column = "control",
                neg.ctl.term = "Neg",
                out.dir = "./output")
# reduce asv table to reasonable asvs
asv_tab_ITS <- asv_tab_ITS[,colSums(asv_tab_ITS) > 99]


# ASSIGN TAXONOMY ####

# for 16S
# if already exists, skip...
if(!file.exists("./output/16S_taxonomy.RDS")){

  tax_16S <- 
    assign_taxonomy_to_asv_table(asv.table = asv_tab_16S,
                                 min.boot = 80,
                                 tax.database = "./taxonomy/silva_nr99_v138.1_wSpecies_train_set.fa.gz")
  saveRDS(tax_16S,"./output/16S_taxonomy.RDS")
} else {
  tax_16S <- readRDS("./output/16S_taxonomy.RDS")
}



# for ITS
# if already exists, skip...
if(!file.exists("./output/ITS_taxonomy.RDS")){
  tax_ITS <- 
    assign_taxonomy_to_asv_table(asv.table = asv_tab_ITS,
                                 min.boot = 80,
                                 tax.database = "./taxonomy/Eukaryome_General_Longread_v1.9_reformatted.fasta.gz")
  saveRDS(tax_ITS,"./output/ITS_taxonomy.RDS")
} else {
  tax_ITS <- readRDS("./output/ITS_taxonomy.RDS")
}


# BUILD PHYSEQ ####
# not including phylogenetic treee yet

# object preparation
# subset metadata
meta_16S <- meta[meta$region == "16S" & (meta$control != "Neg" | is.na(meta$control)),]
meta_ITS <- meta[meta$region == "ITS" & (meta$control != "Neg" | is.na(meta$control)),]

# metadata
met_16S <- sample_data(meta_16S)
met_ITS <- sample_data(meta_ITS)
sample_names(met_16S) <- meta_16S$sample_name
sample_names(met_ITS) <- meta_ITS$sample_name

# ASV tables
otu_16S <- otu_table(asv_tab_16S,taxa_are_rows = FALSE)
otu_ITS <- otu_table(asv_tab_ITS,taxa_are_rows = FALSE)
sample_names(otu_16S) <- meta_16S$sample_name
sample_names(otu_ITS) <- meta_ITS$sample_name

# Tax tables
tax_16S <- tax_table(tax_16S)
tax_ITS <- tax_table(tax_ITS)

# combine parts
ps_16S <- 
  phyloseq(met_16S,otu_16S,tax_16S)
ps_ITS <- 
  phyloseq(met_ITS,otu_ITS,tax_ITS)

# export objects
saveRDS(ps_16S,"./output/16S_Physeq_not-cleaned.RDS")
saveRDS(ps_ITS,"./output/ITS_Physeq_not-cleaned.RDS")

  

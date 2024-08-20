# TRIM ADAPTORS AND RUN QC FILTRATION

# SETUP ####

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

## external software paths (hard coded)
# cutadapt <- "/uufs/chpc.utah.edu/sys/installdir/cutadapt/3.5/bin/cutadapt"

## variables ####
# 16S
bact_27F <- "GAGAGTTTGATCCTGGCTCAG"
bact_1541R <- "AAGGAGGTGATCCAGCCGCA"

# ITS1
ITS9MUNngs <- "TACACACCGCCCGTCG" 
ITS4ngsUni <- "CCTSCSCTTANTDATATGC" 


## load metadata ####
meta <- read_csv("./data/Metadata.csv") %>% 
  janitor::clean_names()

## data checks ####
if(any(!file.exists(meta$filepaths))){
  warning("Some files not present in given paths. Check download. Subsetting to only present files.")
  meta <- meta[file.exists(meta$filepaths),]
}  


# 
# ## setup file paths for QC data
# meta$filtN_paths <- 
#   meta$filepaths %>% 
#   str_replace("data/16S/","data/16S/filtN/") %>% 
#   str_replace("data/ITS/","data/ITS/filtN/") %>% 
#   str_replace(".fastq.gz","_filtN.fastq.gz")
# 
# meta$cutadapt_paths <- 
#   meta$filepaths %>% 
#   str_replace("data/16S/","data/16S/cutadapt/") %>% 
#   str_replace("data/ITS/","data/ITS/cutadapt/") %>% 
#   str_replace(".fastq.gz","_cutadapt.fastq.gz")
# 
# meta$filtered_paths <- 
#   meta$filepaths %>% 
#   str_replace("data/16S/","data/16S/filtered/") %>% 
#   str_replace("data/ITS/","data/ITS/filtered/") %>% 
#   str_replace(".fastq.gz","_filtered.fastq.gz")

# make directories, if needed       
new.dirs <- c("data/16S/filtN/","data/ITS/filtN/",
              "data/16S/cutadapt/","data/ITS/cutadapt/",
              "data/16S/filtered/","data/ITS/filtered/")
for(i in new.dirs){if(!dir.exists(i)){dir.create(i)}}


# PREFILTER ####

# if already done, automatically skip
filtN.fwds <- meta$filepaths[!file.exists(meta$filt_n_paths)]
filtN.outs <- meta$filt_n_paths[!file.exists(meta$filt_n_paths)]

if(length(filtN.fwds) > 0){
  filterAndTrim(fwd = filtN.fwds,
                filt = filtN.outs,
                maxN = 0,
                compress = TRUE,
                multithread = parallel::detectCores())
} else {
  message("All files previously had Ns removed...skipping filtN step.")
}



# RUN CUTADAPT ####

# if already done, automatically skip
if(any(!file.exists(meta$cutadapt_paths))){
  # on 16S
  remove_primers(metadata = meta, 
                 amplicon.colname = "region",
                 amplicon = "16S", 
                 sampleid.colname = "sample_name", 
                 filtN.colname = "filt_n_paths", 
                 cutadapt.colname = "cutadapt_paths",
                 fwd_primer = bact_27F,
                 rev_primer = bact_1541R,
                 multithread = parallel::detectCores()-1)
  
  
  # on ITS
  remove_primers(metadata = meta, 
                 amplicon.colname = "region",
                 amplicon = "ITS", 
                 sampleid.colname = "sample_name", 
                 filtN.colname = "filt_n_paths", 
                 cutadapt.colname = "cutadapt_paths",
                 fwd_primer = ITS9MUNngs,
                 rev_primer = ITS4ngsUni,
                 multithread = parallel::detectCores()-1)
} else {
  message("All cutadapt files already exist... Skipping primer removal step.")
}

# FILTER AND TRIM ####

if(any(!file.exists(meta$filtered_paths))){
  # for 16S data
  qc_filter(metadata = meta,
            amplicon.colname = "region",
            amplicon = "16S",
            input.colname = "cutadapt_paths",
            output.colname = "filtered_paths",
            max.ee = 2,
            trucq = 2)
  
  # for ITS data
  qc_filter(metadata = meta,
            amplicon.colname = "region",
            amplicon = "ITS",
            input.colname = "cutadapt_paths",
            output.colname = "filtered_paths",
            max.ee = 2,
            trucq = 2)
} else {
  message("All QC filtered files already exist... Skipping filter-and-trim step.")
}



  
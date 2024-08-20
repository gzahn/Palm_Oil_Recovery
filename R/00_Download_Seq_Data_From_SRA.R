#------------------------------------------------------------------#
# Download raw sequence data from SRA Project PRJNA1117193         #
#                                                                  #
# These are single-end PacBio reads                                #
# Requires that sra-toolkit be installed and in your system path   #
# (https://github.com/ncbi/sra-tools)                              #
#------------------------------------------------------------------#


# Setup ####
library(tidyverse)
df <- read_csv("./data/Metadata.csv")

# separate out ITS and 16S samples
fung <- df %>% 
  dplyr::filter(region == "ITS")
bact <- df %>% 
  dplyr::filter(region == "16S")

# Download fungal amplicons ####
accessions <- fung$accession
filenames_f <- file.path(str_remove(fung$filename,".gz"))
dl_filenames_1 <- paste0(accessions,"_1.fastq")

# run download in a for-loop for all accessions
# compress and rename and move files
if(!dir.exists("./data/ITS")){
  dir.create("./data/ITS",recursive = TRUE)
}

for(i in seq_along(accessions)){
  system2("fasterq-dump",args = accessions[i])
  system2("gzip", args = c(dl_filenames_1[i]))
  file.rename(paste0(dl_filenames_1[i],".gz"),paste0("./data/ITS/",fung$filename[i]))
}



# Download bacterial amplicons ####
accessions <- bact$accession
filenames_f <- file.path(str_remove(bact$filename,".gz"))
dl_filenames_1 <- paste0(accessions,"_1.fastq")

# run download in a for-loop for all accessions
# compress and rename and move files
if(!dir.exists("./data/16S")){
  dir.create("./data/16S",recursive = TRUE)
}

for(i in seq_along(accessions)){
  system2("fasterq-dump",args = accessions[i])
  system2("gzip", args = c(dl_filenames_1[i]))
  file.rename(paste0(dl_filenames_1[i],".gz"),paste0("./data/ITS/",bact$filename[i]))
}


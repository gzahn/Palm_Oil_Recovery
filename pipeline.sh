#!/bin/bash

# Run this script to reproduce the entire bioinformatics workflow


# download and prepare taxonomy databases
cd ./taxonomy
./build_taxonomy_database.sh
cd -


# download raw sequence files
Rscript ./R/00_Download_Seq_Data_From_SRA.R

# QC of raw sequences
Rscript ./R/01_Trim_Adaptors.R

# Build phyloseq objects
# incl ASV tables, metadata, taxonomic assignments, phytree for 16S
Rscript 02_Build_Physeq_Objects.R




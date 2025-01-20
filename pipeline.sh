#!/bin/bash

# Run this script to reproduce the entire bioinformatics workflow


# download and unzip phylogeny files
wget https://zenodo.org/records/14708616/files/2024.09.phylogeny.id.nwk.gz?download=1 -O phylogeny_info/2024.09.phylogeny.id.nwk.gz
wget https://zenodo.org/records/14708616/files/2024.09.taxonomy.id_w_order-column.tsv.gz?download=1 -O phylogeny_info/2024.09.taxonomy.id_w_order-column.tsv.gz

gunzip phylogeny_info/2024.09.phylogeny.id.nwk.gz
gunzip phylogeny_info/2024.09.taxonomy.id_w_order-column.tsv.gz

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




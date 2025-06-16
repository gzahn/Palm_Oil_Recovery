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
Rscript ./R/01_Trim_Adaptors_and_Filter.R

# build phyloseq objects
# incl ASV tables, metadata, taxonomic assignments, phytree for 16S
Rscript ./R/02_Build_Physeq_Objects.R

# clean up, remove non-target taxa, etc.
Rscript ./R/03_Clean_Physeq_Objects.R

# build phylogenies for 16S and full-length ITS
Rscript ./R/04_Build_Phylogenies.R

# alpha-diversity analyses
Rscript ./R/05_Alpha_Diversity.R

# beta-diversity analyses, GDM, etc.
Rscript ./R/05_Beta_Diversity.R

# build co-ocurrence networks 
Rscript ./R/06_Network_Analyses.R

# find "important" taxa; corncob, random forest, hub taxa
Rscript ./R/07_Important_Taxa.R

# investigate community assembly trends
Rscript ./R/08_Community_Assembly.R

# look at traits of taxa
Rscript ./R/09_Taxa_Traits.R

# build figs of phylogeny over time
Rscript ./R/10_Phylogenetic_Heatmap_figures.R

# conduct FAVA analysis of community compositional variability
Rscript ./R/11_FAVA_Analysis.R



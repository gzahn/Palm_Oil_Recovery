#phylogenetic tree heatmap figures

#SETUP ####

#packages####
library(ape)
library(ggtree)
library(phytools)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggnewscale)
####read in function used to collapse down identical tips that are monophyletic####
#function used to collapse identical tips that are monophyletic
collapse_identical_tips <- function(phy,tip_label){
  matching_tips <- which(phy$tip.label==tip_label)
  # number of tips in tree
  nt <- length(phy$tip.label)
  # Number of tips matching the label
  nm <- length(matching_tips)
  keep <- numeric(nm)
  
  cur_tip <- 1
  while(cur_tip<=nm){
    if(cur_tip == nm){
      keep[cur_tip] <- 1
      break
    }
    next_tip <- cur_tip + 1
    mrca_ <- getMRCA(phy,c(matching_tips[cur_tip],matching_tips[next_tip]))
    descendants <- getDescendants(phy, mrca_)
    descendant_tips <- descendants[descendants<=nt]
    if(all(descendant_tips %in% matching_tips)){
      keep[cur_tip] <- 1
      cur_tip <- cur_tip + length(descendant_tips)
    }else{
      keep[cur_tip] <- 1
      cur_tip <- cur_tip + 1
    }
  }
  to_drop <- matching_tips[!keep]
  new_phy <- drop.tip(phy,to_drop)
  return(new_phy)
}
####read in bact data####
bact<-read.csv("./output/Bacterial_order_presence_over_time.csv")

#read file in from website
#this file was made by downloading the 2024.09.taxonomy.id.tsv.gz file from
#https://ftp.microbio.me/greengenes_release/current/
#the file was then modified to make a column with just the order information
#the bash script that was used was:
#cat 2024.09.taxonomy.id.tsv | cut -d ";" -f 4 | sed 's/o__$/NA/' | sed 's/o__//' | sed 's/^ //'
#this made a file called 2024.09.taxonomy.id_w_order-column.tsv
bact.names<-read.delim("./phylogeny_info/2024.09.taxonomy.id_w_order-column.tsv", header=T)

#read tree in from website
#this was downloaded from https://ftp.microbio.me/greengenes_release/current/
#this file is the 2024.09.phylogeny.id.nwk file
bact.tree<-read.tree("./phylogeny_info/2024.09.phylogeny.id.nwk")
#check if rooted
#is.rooted(bact.tree)

####changing names on the tips of the bacteria phylogeny
#setting the name map info so we can make changes to the tip labels to reflect
#orders
name_map <- setNames(bact.names$Order, bact.names$Feature.ID)

#changing the names on the tree
bact.tree$tip.label <- name_map[bact.tree$tip.label]
na_tips <- which(is.na(bact.tree$tip.label))

#we want to get rid of tips which are NA (do not have order information)
#something like 2.3 million tips
tips_to_remove <- bact.tree$tip.label[is.na(bact.tree$tip.label) | bact.tree$tip.label == "NA"]
cleaned_tree.bact <- drop.tip(bact.tree, tips_to_remove)
#swap out names for orders which are labelled differently in the tree than we have in the data
#this usually means a letter or a number is added after the order name
#these were identified doing the below and then manually searching the names for 
#the below transformations
#
# setdiff(unique(bact$Order), unique(ult.tree$tip.label))
# setdiff(unique(bact.names$Order), unique(bact$Order))
# toomanynames<-setdiff(unique(bact.names$Order), unique(bact$Order))
# write.csv(toomanynames, "overfloworders.csv")
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Azospirillales_507929"]<-"Azospirillales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Bacillales_K"]<-"Bacillales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Cytophagales_B"]<-"Cytophagales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Enterobacterales_639860"]<-"Enterobacterales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Flavobacteriales_B_894219"]<-"Flavobacteriales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Haliangiales_463188"]<-"Haliangiales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Holosporales_A"]<-"Holosporales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Pseudomonadales_A_628887"]<-"Pseudomonadales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Rhizobiales_505101"]<-"Rhizobiales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Rhodospirillales_A_505895"]<-"Rhodospirillales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Rickettsiales_A"]<-"Rickettsiales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Spirochaetales_E"]<-"Spirochaetales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Streptomycetales_400645"]<-"Streptomycetales"
cleaned_tree.bact$tip.label[cleaned_tree.bact$tip.label =="Thermoactinomycetales_367867"]<-"Thermoactinomycetales"


####additional trimming of the bact tree to fit our dataset####
#remove orders in tree but not in dataset
#this trims down to about 12.9 million tips
only.in.tree<- setdiff(unique(cleaned_tree.bact$tip.label), unique(bact$Order))
cleaned_tree.bact<-drop.tip(cleaned_tree.bact, tip=only.in.tree)

####cleaning up the bact tree####
#make a list of duplicated tip names
duplicated_names <- unique(cleaned_tree.bact$tip.label[duplicated(cleaned_tree.bact$tip.label)])
#duplicated_names

#i know this will collapse down to the correct number of tips, so I use this to
#do it all at once
while(length(cleaned_tree.bact$tip.label)>length(unique(cleaned_tree.bact$tip.label))){
 for(i in 1:length(duplicated_names)){
   cleaned_tree.bact<-collapse_identical_tips(cleaned_tree.bact, duplicated_names[i])
print(i)
 }
 }
saveRDS(cleaned_tree.bact,"./phylogeny_info/cleaned_tree_bact.RDS")
#you can also use the below until you find a tip number you want.  It runs
#six times
# for(i in 1:length(duplicated_names)){
#        cleaned_tree.bact<-collapse_identical_tips(cleaned_tree.bact, duplicated_names[i])
#        print(i)
#      }
#the below will plot the tree on its own
#plot(cleaned_tree.bact,cex=0.7)

#make tree ultrametric
ult.tree.bact<-chronos(cleaned_tree.bact)
#here's the tree
plot(ult.tree.bact)

####working with the bact tree####
#tree name is ult.tree.bact
#data is called bact.dat

#filter to data that is availalbe for the tree
filtered_bact <- bact %>% filter(Order %in% ult.tree.bact$tip.label)

#filter to only the data we want
NASV.bact <- filtered_bact %>% select(Order, time, N_ASVs)
#make year a factor and fix its labels
NASV.bact$time<-as.factor(NASV.bact$time)
levels(NASV.bact$time)<-c("Oil palm", "0-1 Years", "1-2 Years", "2-3 Years", "3-4 Years")
#fix the formating
bact.abundance_data <- NASV.bact %>%
  pivot_wider(names_from = time, values_from = N_ASVs)

#this makes it so we can have it as a matrix for phytools
bact.abundance_data<-as.data.frame(bact.abundance_data)
rownames(bact.abundance_data) <- bact.abundance_data$Order
bact.abundance_data <- bact.abundance_data[, -1]

#this is for relative abundance in bacteria
rel.bac <- filtered_bact %>% select(Order, time, relative_abundance)
#make year a function and fix its labels
rel.bac$time<-as.factor(rel.bac$time)
levels(rel.bac$time)<-c("Oil palm", "0-1 Years", "1-2 Years", "2-3 Years", "3-4 Years")
#fix the formating
rel_abundance_data.bact <- rel.bac %>%
  pivot_wider(names_from = time, values_from = relative_abundance)
#this makes it so we can have it as a matrix for phytools
rel_abundance_data.bact<-as.data.frame(rel_abundance_data.bact)
rownames(rel_abundance_data.bact) <- rel_abundance_data.bact$Order
rel_abundance_data.bact <- rel_abundance_data.bact[, -1]

abund.mat.bact<-as.matrix(bact.abundance_data)
rel.abund.mat.bact<-as.matrix(rel_abundance_data.bact)


#plots the tree
p <- ggtree(ult.tree.bact) + 
  geom_tiplab(size=2.75, align=TRUE) + 
  theme(plot.margin = margin(t = 8, r = 0.2, b = 15, l = 1)) 
#first heatmap of number of unique ASVs
p1 <- gheatmap(p, abund.mat.bact, offset=.4, width=.3,
               colnames_angle=60, colnames_offset_y =-3,colnames_offset_x = -.03,
               color="black") +
  scale_fill_viridis_c(option="D", direction=-1,
                       end=0.9,name="# Unique\nASVs",
                       na.value = "white")
p2 <- p1 + new_scale_fill()
#second heatmap
bact_final_plot <- 
  gheatmap(p2, rel.abund.mat.bact, offset=0.75, width=.3,
         colnames_angle=60, colnames_offset_y = -3,colnames_offset_x = -.03,
         colnames_position = "bottom",
         color="black") +
  scale_fill_viridis_c(option="A", direction=-1,
                       end=0.9, name="Relative\nAbundance",
                       na.value = "white")+
  #the ylim max is important for the number of orders we have
  ylim(NA, max(ult.tree.bact$edge.length) + 100)+
  annotate(geom="text", x=1.43, y=100, label="ASVs",
           color="black", hjust=0)+
  annotate(geom="text", x=1.78, y=100, label="Rel. Abund.",
           color="black", hjust=0)
bact_final_plot
ggsave("./output/figs/bact_tree_heatmap.png",dpi=400, height = 18,width = 10)

#Fungus tree stuff####
#if you're only working with this stuff, be sure to read in the function at the top of the script

#read in tree and data####
#read in tree from paper
#from this figshare:https://figshare.com/articles/dataset/Scripts_and_analyses_used_for_the_fungal_phylogeny/12751736
#the file is 1672taxa_290genes_bb_1
fungi.tree<-read.tree("./phylogeny_info/1672taxa_290genes_bb_1.treefile")
#read in fungus data
fung.dat<-read.csv("./output/Fungal_order_presence_over_time.csv")
#read in the translation information to go from species tip labels 
###to order names as indicated by supp data in file
fun.namelist<-read.csv("./phylogeny_info/fung.tax.names.csv")

#make modifications to the tree to fit our data
####make a list of species indicated as outgroups to root the tree####
out.group<-c("Fonticula_alba",
             "Fonticula_like_sp_SCN_57_25",
             "Capsaspora_owczarzaki_ATCC_30864",
             "Corallochytrium_limacisporum",
             "Creolimax_fragrantissima",
             "Ichthyophonus_hoferi",
             "Ichthyosporea_sp._XGB_2017a",
             "Monosiga_brevicollis_MX1",
             "Salpingoeca_rosetta",
             "Sphaeroforma_arctica",
             "Sphaeroforma_sirkka",
             "Strigamia_maritima",
             "Mnemiopsis_leidyi",
             "Helobdella_robusta",
             "Lingula_anatina",
             "Ixodes_scapularis",
             "Branchiostoma_lanceolatum",
             "Capitella_teleta",
             "Trichoplax_adhaerens",
             "Ciona_intestinalis",
             "Daphnia_pulex",
             "Nematostella_vectensis",
             "Octopus_bimaculoides",
             "Lottia_gigantea",
             "Mus_musculus",
             "Drosophila_melanogaster",
             "Strongylocentrotus_purpuratus",
             "Amphimedon_queenslandica")

#####root the tree####
fungi.tree.root<-root(fungi.tree, outgroup=out.group, resolve.root = T)
#setting thename map info
fun.name_map <- setNames(fun.namelist$Order_linked_to_NCBI, fun.namelist$new_TaxonID)
#changing the names on the tree
fungi.tree.root$tip.label <- fun.name_map[fungi.tree.root$tip.label]

####cleaning up the tree####
#remove outgroups
fun.cleaned_tree<-drop.tip(fungi.tree.root, tip=out.group)
#remove "no ranks:
fun.cleaned_tree<-drop.tip(fun.cleaned_tree, tip="no_rank")
fun.na_tips <- which(is.na(fun.cleaned_tree$tip.label))
fun.cleaned_tree <- drop.tip(fun.cleaned_tree, fun.cleaned_tree$tip.label[fun.na_tips])

#remove orders in tree but not in dataset
fun.only.in.tree<- setdiff(unique(fun.cleaned_tree$tip.label), unique(fung.dat$Order))
fun.cleaned_tree<-drop.tip(fun.cleaned_tree, tip=fun.only.in.tree)

#make a list of duplicated tip names
fun.duplicated_names <- unique(fun.cleaned_tree$tip.label[duplicated(fun.cleaned_tree$tip.label)])
#duplicated_names

#collapsing down the tree's extra tips
while(length(fun.cleaned_tree$tip.label)>length(unique(fun.cleaned_tree$tip.label))){
  for(i in 1:length(fun.duplicated_names)){
    fun.cleaned_tree<-collapse_identical_tips(fun.cleaned_tree, fun.duplicated_names[i])  
  print(i)
    }
}

#make ultrametric
fun.ult.tree<-chronos(fun.cleaned_tree)
#here's the tree
plot(fun.ult.tree)

####working with the tree####
#tree name is fun.ult.tree
#data is called fung.dat

#filter to data that is available for the tree
filtered_fung <- fung.dat %>% filter(Order %in% fun.ult.tree$tip.label)

#filter to only the data we want
NASV.fun <- filtered_fung %>% select(Order, time, N_ASVs)
#make year a factor and fix its labels
NASV.fun$time<-as.factor(NASV.fun$time)
levels(NASV.fun$time)<-c("Oil palm", "0-1 Years", "1-2 Years", "2-3 Years", "3-4 Years")
#fix the formating
fun.abundance_data <- NASV.fun %>%
  pivot_wider(names_from = time, values_from = N_ASVs)
#this makes it so we can have it as a matrix for phytools
fun.abundance_data<-as.data.frame(fun.abundance_data)
rownames(fun.abundance_data) <- fun.abundance_data$Order
fun.abundance_data <- fun.abundance_data[, -1]

#this is for relative abundance
rel.fun <- filtered_fung %>% select(Order, time, relative_abundance)
#make year a function and fix its labels
rel.fun$time<-as.factor(rel.fun$time)
levels(rel.fun$time)<-c("Oil palm", "0-1 Years", "1-2 Years", "2-3 Years", "3-4 Years")
#fix the formating
fun.rel_abundance_data <- rel.fun %>%
  pivot_wider(names_from = time, values_from = relative_abundance)
#this makes it so we can have it as a matrix for phytools
fun.rel_abundance_data<-as.data.frame(fun.rel_abundance_data)
rownames(fun.rel_abundance_data) <- fun.rel_abundance_data$Order
fun.rel_abundance_data <- fun.rel_abundance_data[, -1]
# 
# 

fun.abund.mat<-as.matrix(fun.abundance_data)
fun.rel.abund.mat<-as.matrix(fun.rel_abundance_data)


#makes the tree
f.p <- ggtree(fun.ult.tree) + 
  geom_tiplab(size=3, align=TRUE) + 
  theme(plot.margin = margin(t = 25, r = 0.2, b = 50, l = 1)) 
#makes the first heatmap for # of unique ASVs
f.p1 <- gheatmap(f.p, fun.abund.mat, offset=.25, width=.3,
               colnames_angle=60, colnames_offset_y =-3,colnames_offset_x = -.03,
               color="black") +
  scale_fill_viridis_c(option="D", direction=-1,
                       end=0.9,name="# Unique\nASVs",
                       na.value = "white")
f.p2 <- f.p1 + new_scale_fill()
#second heatmap
fung_final_plot <- 
  gheatmap(f.p2, fun.rel.abund.mat, offset=0.6, width=.3,
         colnames_angle=60, colnames_offset_y = -3,colnames_offset_x = -.03,
         colnames_position = "bottom",
         color="black") +
  scale_fill_viridis_c(option="A", direction=-1,
                       end=0.9, name="Relative\nAbundance",
                       na.value = "white")+
  ylim(NA, max(fun.ult.tree$edge.length) + 40)+
  annotate(geom="text", x=1.28, y=37, label="ASVs",
           color="black", hjust=0)+
  annotate(geom="text", x=1.63, y=37, label="Rel. Abund.",
           color="black", hjust=0)
fung_final_plot
ggsave("./output/figs/fung_tree_heatmap.png",dpi=400,height = 18,width = 10)

library(patchwork)
bact_final_plot | fung_final_plot
ggsave("./output/figs/combined_phylogeny_heatmap.png",dpi=400,height = 18,width = 20)

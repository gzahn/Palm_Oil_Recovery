# palettes ####
pal <- list(
  pal.earthtones = c("#4E6172","#D57500","#8F3B1B","#404F24","#613318","#668D3C"),
  pal.okabe = c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
)
# trt.pal <- c("#b80006","#b36d1e","#e0cb2d","#8ecc12","#27c4a5")
trt.pal <- c("#300206","#b80006","#d65002","#929412","#1A9850")


# %ni% ####
# "not in" logical expression
'%ni%' <- Negate('%in%')


# ra() ####
# relative abundance transformation
ra <- function(x){x/sum(x)}

# pa() ####
# presence/absence transformation
pa <- function(x){ifelse(x>0,1,0)}


# EE() ####
# calculate expected errors in a sequence from a vector of quality scores
EE <- function(qual.scores){
  sum(10^(-qual.scores/10))
}

# plot_bar2() ####
# phyloseq bar plot without lines for each ASV
plot_bar2 <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
          title = NULL, facet_grid = NULL, width = 0.9) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack", width = width)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

# remove_primers() ####
# Function to remove primers from raw amplicon files

remove_primers <- function(metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
                           amplicon.colname = "region", # column name that contains the amplicon info for each sample
                           amplicon = "16S", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
                           sampleid.colname = "sample_name", # column name in metadata containing unique sample identifier
                           filtN.colname = "filt_n_paths", # name of column in metadata indicating fwd filepath to raw data
                           cutadapt.colname = "cutadapt_paths", # name of column with output filenames for cutadapt
                           fwd_primer="GTGCCAGCMGCCGCGGTAA",
                           rev_primer="GGACTACHVGGGTWTCTAAT",
                           multithread=parallel::detectCores()-1){

  library(tidyverse); packageVersion("tidyverse")
  library(dada2); packageVersion("dada2")
  library(Biostrings,lib.loc = mylib); packageVersion("Biostrings")
  library(ShortRead,lib.loc = mylib); packageVersion("ShortRead")
  library(parallel); packageVersion("parallel")
  
  
  # File parsing
  
  # subset metadata to just that amplicon and get file paths
  x <- metadata[metadata[[amplicon.colname]] == amplicon,]
  
  # pull out file paths for convenience
  
  fnFs.filtN <- x[[filtN.colname]]
  fnFs.cutadapt <- x[[cutadapt.colname]]
  # get sample_names
  sample_names <- x[[sampleid.colname]]
  
  FWD <- fwd_primer # Sequence of FWD primer
  REV <- rev_primer  # Sequence of REV primer
  
  allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
  }
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  
  # get revcomplements of primers
  FWD.RC <- dada2:::rc(FWD)
  REV.RC <- dada2:::rc(REV)
  
  
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  R1.flags <- paste("-g", FWD, "-a", REV.RC) 
  # Run Cutadapt
  for(i in seq_along(fnFs.filtN)) {
    print(fnFs.filtN[i])
      system2("cutadapt", args = c(R1.flags,
                                   "-n", 2, # -n 2 required to remove FWD and REV from reads
                                   "--minimum-length 500", # for long reads 
                                   "--cores 0",
                                   "-o", fnFs.cutadapt[i], # output files
                                   fnFs.filtN[i])) # input files
    }
  }

# qc_filter() ####
# run filter_and_trim
qc_filter <- 
  function(metadata,
           amplicon.colname = "region",
           amplicon = "16S",
           input.colname = "cutadapt_paths",
           output.colname = "filtered_paths",
           max.ee = 1,
           trucq = 2,
           threads = parallel::detectCores()-1){
    
    x <- metadata[metadata[[amplicon.colname]] == amplicon,]
    
    # pull out file paths for convenience
    input.files <- x[[input.colname]]
    output.files <- x[[output.colname]]
    
    out <- filterAndTrim(fwd = input.files, 
                         filt = output.files,
                         maxN=0,
                         maxEE=max.ee,
                         truncQ=trucq,
                         rm.phix=TRUE,
                         compress=TRUE,
                         multithread=threads)
    saveRDS(out,paste0("./data/",amplicon,"/filtration_stats.RDS"))
  }




# build_asv_table() ####
# function to run dada2 pipeline on trimmed amplicon reads
# error profiling and correction should be done by sequencing run
# this function will run dada2 on files from a single sequencing run, given a metadata sheet that has samples from multiple runs

build_asv_table <- 
  function(metadata = meta,
           filtered_path = "./data/16S/filtered",
           amplicon = "16S",
           multithread = parallel::detectCores()-1,
           control.column = "control",
           neg.ctl.term = "Neg",
           out.dir = "./output"){
    
    if(!dir.exists(out.dir)){
      dir.create(out.dir,recursive = TRUE)
    }
    
    
    
    # get existing filtered files
    filts_f <- list.files(filtered_path,full.names = TRUE,pattern = "fastq.gz")
    
    # learn errors
    
    # if already done, skip...
    errF_out <- paste0("ErrorModel","_",amplicon,".RDS")
    if(!file.exists(file.path(out.dir,errF_out))){
      errF <- learnErrors(fls = filts_f, 
                          multithread = ifelse(multithread>1,TRUE,FALSE), 
                          MAX_CONSIST = 20,
                          verbose = TRUE,
                          randomize = TRUE) # set multithread = FALSE on Windows
      
      saveRDS(errF,file.path(out.dir,errF_out))
    } else {
      errF <- readRDS(file.path(out.dir,errF_out))
    }
    
    
       
    # if already done, skip...
    if(!file.exists(file.path(out.dir,paste0("DADA_",amplicon,".RDS")))){
      # SAMPLE INFERRENCE 
      
      # run in for-loop to reduce memory requirements
      for(sam in seq_along(filts_f)){
        # dereplicate seqs
        derep <- derepFastq(filts_f[sam], verbose=TRUE)
        cat("Processing: ", filts_f[sam], "\n")
        dadaFs[sam] <- dada(derep=derep, 
                            err=errF, 
                            multithread=ifelse(multithread>1,TRUE,FALSE),
                            verbose=TRUE) 
      }
      saveRDS(dadaFs,file.path(out.dir,paste0("DADA_",amplicon,".RDS")))
    } else {
      dadaFs <- readRDS(file.path(out.dir,paste0("DADA_",amplicon,".RDS")))
    }
    
    
    # Make seq table 
    seqtab <- makeSequenceTable(dadaFs)
    
    # if already done, skip...
    if(!file.exists(paste0(out.dir,"/",amplicon,"_ASV_Table.RDS"))){
      # remove chimeras
      seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                          method="consensus", 
                                          multithread=ifelse(multithread>1,TRUE,FALSE), 
                                          verbose=TRUE)
      
      # track read counts so far
      # getN <- function(x) sum(getUniques(x))
      # # find read counts for raw and filtered data files
      # meta$filepaths
      # track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
      # colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
      # rownames(track) <- sample_names
      # # make name for read tracking output
      # track.out <- paste0(out.dir,"/",amplicon,"_trackreads.RDS")
      # # export tracking results
      # saveRDS(track,track.out)
      
      # REMOVE CONTAMINANTS 
      
      # find negative control samples, if any
      metadata[[control.column]][is.na(metadata[[control.column]])] <- FALSE
      metadata[["neg_ctl"]] <- ifelse(metadata[[control.column]] == neg.ctl.term,
                                      TRUE,FALSE)
      # subset metadata to just this amplicon
      metadata <- metadata[metadata[['region']] == amplicon,]
      
      # only run if there are negative control(s) that have at least some reads
      if(any(metadata[["neg_ctl"]]) & sum(seqtab.nochim[metadata[["neg_ctl"]],]) > 0){
        # Find and remove contaminants
        contams = decontam::isContaminant(seqtab.nochim, neg = metadata[["neg_ctl"]], normalize = TRUE)
        
        # remove contaminant sequences and control samples from both tables, respectively ####
        seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
        seqtab.nochim = seqtab.nochim[!metadata[["neg_ctl"]],]
        print(paste0(sum(contams$contaminant)," likely contaminants were removed from the ASV table of ",amplicon,"."))
      }
      
      # make output name for ASV table
      asv_out <- paste0(out.dir,"/",amplicon,"_ASV_Table.RDS")
      # export
      saveRDS(seqtab.nochim,asv_out)
      
    } else {
      seqtab.nochim <- readRDS(paste0(out.dir,"/",amplicon,"_ASV_Table.RDS"))
    }
    
    return(seqtab.nochim)
  }


# assign_taxonomy_to_asv_table() ####
# function to assign taxonomy to an asv table from the dada2 pipeline
# inputs: asv table | path-to-database | 

assign_taxonomy_to_asv_table <- function(asv.table, # asv table object name
                                         tax.database, # path to taxonomic database (fasta)
                                         multithread=(parallel::detectCores()-1), # set to FALSE on Windows
                                         random.seed=666,
                                         try.rc = TRUE, # attempt revComplement assignments as well? (doubles time)
                                         min.boot=50 # bootstrap of 50% recommended for seqs shorter than 250nt
                                         ){
  x <- assignTaxonomy(seqs = asv.table,
                      refFasta = tax.database,
                      minBoot = min.boot,
                      tryRC = try.rc,
                      outputBootstraps = FALSE,
                      multithread = multithread,
                      verbose = FALSE)
  
  return(x)
  
}


# find_gps_dists() ####
# given two data.frames (x,y) of lat/lon points, finds the distance from each point 
# in x to each point in y. By default, returns the distance to only the closest point in y
# good for finding nearest points
find_gps_dists <- 
  function(points1,points2,min.only=TRUE){
    
    # tests
    stopifnot(any(class(points1) == "data.frame"))
    stopifnot(any(class(points2) == "data.frame"))
    
    if(ncol(points1) != 2 | ncol(points2) != 2){
      stop("data frames must have 2 columns only")
    }
    
    if(!apply(points1,2,class) %>% unique() %in% c("numeric","integer")){
      stop("columns must be numeric; col1=longitude,col2=latitude")
    }
    
    if(!apply(points2,2,class) %>% unique() %in% c("numeric","integer")){
      stop("columns must be numeric; col1=longitude,col2=latitude")
    }
    
    # actual function
    mylist <- list()
    for(i in 1:nrow(points1)){
      mylist[[i]] <- points1[i,] %>% unlist
    }
    
    distances <- list()
    for(i in 1:nrow(points2)){
      
      mydistfunction <- function(x){geosphere::distHaversine(x,points2[i,])}
      colname <- paste0("dist_to_",i)
      distances[[colname]] <- map_dbl(mylist, mydistfunction)
    }
    
    x <- as.data.frame(distances)
    mins <- apply(x,1,min)
    
    if(min.only){
      return(mins)
    } else {
      return(x)
    }
    
  }


# auto_permanova() ####
auto_permanova <- function(physeq,
                            data,
                            pred.cols,
                            strata.col,
                            mod.type = "additive"){
  
  ra_comm <- 
    physeq %>% 
    transform_sample_counts(ra) %>% 
    otu_table() %>% 
    as.matrix()
  
  cc <- data %>% 
    dplyr::select(all_of(c(pred.cols,strata.col))) %>% 
    complete.cases()
  
  df <- data[cc,]

  if(mod.type == "additive"){
    mod.formula <- as.formula(paste0("ra_comm[cc,]","~",paste(pred.cols,collapse=" + ")))
  }
  
  if(mod.type == "interactive"){
    mod.formula <- as.formula(paste0("ra_comm[cc,]","~",paste(pred.cols,collapse=" * ")))
  }
  
  
  # run simple permanova
  if(!is.na(strata.col)){
    mod <- 
      adonis2(data = df,
              formula = mod.formula,
              strata = df[[strata.col]])
  } else {
    mod <- 
      adonis2(data = df,
              formula = mod.formula)
  }
  
  
 return(mod) 
}


# clean_ps_taxonomy() ####
# get rid of the annoying "k__" stuff at the beginning of taxonomy assignments for each tax level
# these are usually only an issue with fungal assignments (e.g., UNITE format)


clean_ps_taxonomy <- function(physeq,
                              n.ranks=7 # number of taxonomic ranks in physeq object
                              ){
  
  # get rank names
  ranks <- rank_names(physeq)
  prefix <- str_sub(ranks,end=1) %>% str_to_lower() %>% paste0("__")
  
  x <- tax_table(physeq)@.Data
  
  for(i in seq_along(ranks)){
    x[,i] <- x[,i] %>% str_remove(prefix[i])
  }
  
  out <- phyloseq(otu_table(physeq,taxa_are_rows = FALSE),
                  tax_table(x),
                  sample_data(physeq))
  
  return(out)
}


# simplify_fungal_guilds()
# extract major guild groupings
# needs a data.frame with a "Guild" column that contains the results from FunGuild assignment
# will return the data.frame with a new column "major_guild"
simplify_fungal_guilds <- 
  function(x){
    x %>% 
      mutate(major_guild = case_when(grepl("Ectomycorrhizal",Guild,ignore.case = TRUE) ~ "Ectomycorrhizal",
                                     grepl("ericoid",Guild,ignore.case = TRUE) ~ "Ericoid mycorrhizal",
                                     grepl("arbuscular",Guild,ignore.case=TRUE) ~ "Arbuscular mycorrhizal",
                                     grepl("Plant Pathogen",Guild,ignore.case=TRUE) ~ "Plant pathogen",
                                     grepl("Animal Pathogen",Guild,ignore.case=TRUE) ~ "Animal pathogen",
                                     grepl("Saprotroph",Guild,ignore.case=TRUE) &
                                       !grepl("mycorrhizal",Guild,ignore.case=TRUE) &
                                       !grepl("pathogen",Guild,ignore.case=TRUE) ~ "Saprotroph",
                                     grepl("Orchid Mycorrhizal",Guild) ~ "Orchid Mycorrhizal",
                                     grepl("lichenized",Guild,ignore.case=TRUE) ~ "Lichenized",
                                     grepl("Animal Parasite",Guild) ~ "Animal Parasite",
                                     grepl("Algal Parasite|Plant Parasite",Guild) ~ "Plant Parasite",
                                     grepl("Animal Symbiotroph",Guild) ~ "Animal Symbiotroph"
      ))
  }

# googlemap_json_to_string()
# convert json google map styling to api string
googlemap_json_to_string <- 
  function (style_list) 
  {
    style_string <- ""
    for (i in 1:length(style_list)) {
      if ("featureType" %in% names(style_list[[i]])) {
        style_string <- paste0(style_string, "feature:", 
                               style_list[[i]]$featureType, "|")
      }
      if ("elementType" %in% names(style_list[[i]])) {
        style_string <- paste0(style_string, "element:", 
                               style_list[[i]]$elementType, "|")
      }
      elements <- style_list[[i]]$stylers
      a <- lapply(elements, function(x) paste0(names(x), ":", 
                                               x)) %>% unlist() %>% paste0(collapse = "|")
      style_string <- paste0(style_string, a)
      if (i < length(style_list)) {
        style_string <- paste0(style_string, "&style=")
      }
    }
    style_string <- gsub("#", "0x", style_string)
    style_string <- gsub("[|]", "%7C", style_string)
    return(style_string)
  }


# Function to get the last two words
get_last_two_words <- function(x) {
  words <- unlist(strsplit(x, " "))
  last_two <- tail(words, 2)
  paste(last_two, collapse = " ")
}

# pull out multcomp letters for statistical groupings
# this ONLY works for "treatment" in this Palm Oil data set
multcomp_letters <- 
  function(mod){
    tuk <- TukeyHSD(mod)
    cld <- multcompLetters4(mod,tuk)
    df <- 
      cld %>% 
      pluck('treatment') %>% 
      pluck('Letters') %>% 
      as.data.frame()
    df['treatment'] <- row.names(df)
    names(df)[1] <- 'group'
    df$treatment <- factor(df$treatment,levels=levels(meta$treatment))
    
    return(df)
  }

# treatment vs metadata plot
#ONLY works for this Palm Oil data set
treatment_meta_plot <- 
  function(y){
    ggboxplot(meta, x = "treatment", y = y,
              fill = "treatment", palette = trt.pal,
              add = "jitter") +
      labs(x="\nTreatment") +
      theme(legend.position = 'none',
            line = element_line(linewidth = 6),
            axis.title = element_text(face='bold',size=18),
            axis.text = element_text(face='bold',size=14))
  }


# Get igraph network attributes
find_ig_subset_attr <- function(ps.subset,ig.full){
  
  # run checks and tests for function
  stopifnot(class(ps.subset) == "phyloseq")
  stopifnot(class(ig.full) == "igraph")
  
  # get present taxa in physeq subset
  present_vertices <- which(taxa_sums(ps.subset) > 0)
  
  # build subgraph from present taxa only
  ps.subset %>% ntaxa()
  ig.subset <- igraph::subgraph(graph = ig.full,vids = present_vertices)
  # calculate various ig params
  # mean_alpha_centrality <- mean(igraph::alpha_centrality(ig.subset,tol=0),na.rm = TRUE)
  clique_num <- igraph::clique_num(ig.subset)
  mean_betweenness <- mean(igraph::betweenness(ig.subset),na.rm = TRUE)
  mean_closeness <- mean(igraph::closeness(ig.subset),na.rm = TRUE)
  mean_coreness <- mean(igraph::coreness(ig.subset),na.rm = TRUE)
  deg_dist <- igraph::degree_distribution(ig.subset)
  global_effic <- igraph::global_efficiency(ig.subset)
  n_vertices <- igraph::vcount(ig.subset)
  n_edges <- igraph::ecount(ig.subset)
  mean_dist <- igraph::mean_distance(ig.subset)
  similarity_matrix <- igraph::similarity(ig.subset,method='jaccard')
  clustering_coeficient <- igraph::transitivity(ig.subset)
  max_degree <- max(igraph::degree(ig.subset),na.rm=TRUE)
  
  out <- list(ig=ig.subset,
              n_vertices=n_vertices,
              n_edges=n_edges,
              mean_dist=mean_dist,
              # mean_alpha_centrality=mean_alpha_centrality,
              clique_num=clique_num,
              mean_betweenness=mean_betweenness,
              mean_closeness=mean_closeness,
              mean_coreness=mean_coreness,
              deg_dist=deg_dist,
              global_effic=global_effic,
              similarity_matrix=similarity_matrix,
              clustering_coeficient=clustering_coeficient,
              max_degree=max_degree)
  
  
  return(out)
}

# is_increasing() ####
# determines whether a numeric vector consistently increases in value along its whole length
is_increasing <- function(vec) {
  return(all(diff(vec) > 0))
}

# is_greater_than_first() ####
is_greater_than_first <- function(vec){
  first_element <- vec[1]
  other_elements <- vec[2:length(vec)]
  return(all(other_elements > first_element))
}


# get_bacdive_morphology() ####
# bacdive helper function
get_bacdive_morphology <- function(genus){
  x <- BacDive::request(object = bacdive,
                        query = genus,
                        search = "taxon")
  # find bacdive id
  y <- x$results
  
  if(length(y) == 0){
    df <- data.frame(genus=genus,
                     length=NA,
                     width=NA,
                     shape=NA)
    return(df)
  }
  
  # get info on that id
  z <- BacDive::fetch(bacdive,y)
  
  # get morphology
  m <- z$results %>% 
    map("Morphology")
  
  m <- m %>% map("cell morphology")
  # subset to non-null entries
  m <- m[c(which(m %>% map(notnull) %>% unlist) %>% names)]
  
  # test to see if morphology present
  # if morphology of type strain is empty, check taxonomic synonyms
  if(length(m) == 0){
    
    # find synonyms
    taxonomy <- z$results %>% 
      map("Name and taxonomic classification")
    syn <- taxonomy %>% 
      map("LPSN") %>% 
      map("synonyms") %>% 
      map("synonym") %>% 
      unlist() %>% 
      unique()
    
    # if not synonyms, return NAs
    if(is.null(syn)){
      df <- data.frame(genus=genus,
                       length=NA,
                       width=NA,
                       shape=NA)
      return(df)
    } else {
      # re-run bacdive request with synonyms
      x <- BacDive::request(object = bacdive,
                            query = syn,
                            search = "taxon")  
      
      if(length(x$results) == 0){
        df <- data.frame(genus=genus,
                         length=NA,
                         width=NA,
                         shape=NA)
        return(df)
      } else {
        z <- BacDive::fetch(bacdive,x$results)
      }
      
    }
  }
  
  # get morphology again (perhaps using synonyms)
  m <- z$results %>% 
    map("Morphology")
  m <- m %>% map("cell morphology")
  # subset to non-null entries
  m <- m[c(which(m %>% map(notnull) %>% unlist) %>% names)]
  
  # test to see if morphology present
  # if morphology still not present, export NAs
  
  if(length(m) == 0){
    df <- data.frame(genus=genus,
                     length=NA,
                     width=NA,
                     shape=NA)
    return(df)
  } else {
    # get morphology
    m <- z$results %>% 
      map("Morphology")
    m <- m %>% map("cell morphology")
    # subset to non-null entries
    m <- m[c(which(m %>% map(notnull) %>% unlist) %>% names)]
    
    # cell length
    length <- m %>% 
      map("cell length")
    if(which(length %>% map(notnull) %>% unlist()) %>% length() ==0){
      length <- NA
    } else {
      length <- length[which(length %>% map(notnull) %>% unlist())] %>% 
        unlist()
    }
    
    # cell width
    width <- m %>% 
      map("cell width")
    if(which(width %>% map(notnull) %>% unlist()) %>% length() ==0){
      width <- NA
    } else {
      width <- width[which(width %>% map(notnull) %>% unlist())] %>% 
        unlist()
    }
    
    # cell shape
    shape <- m %>% 
      map("cell shape")
    if(which(shape %>% map(notnull) %>% unlist()) %>% length() ==0){
      shape <- NA
    } else{
      shape <- shape[which(shape %>% map(notnull) %>% unlist())] %>% 
        unlist()
    }
    
    l <- length(length)
    w <- length(width)
    s <- length(shape)
    max_rows <- min(c(l,w,s))
    
    length <- length[1:max_rows]
    width <- width[1:max_rows]
    shape <- shape[1:max_rows]
    
    df <- data.frame(genus=genus,
                     length=length,
                     width=width,
                     shape=shape)
    return(df)
  }
  
}
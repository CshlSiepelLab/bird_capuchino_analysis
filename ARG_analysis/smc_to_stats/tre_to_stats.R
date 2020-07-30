  ### R script to run on Linux  ###
############################################################################
##### read tree file and compute stats for every tree                   ####
####  Stats computed for every tree:                                    ####
####  - clade60age - age of youngest clade of size >= 60                ####
####                 (proxy for general depth of clades in tree)        ####
####  - clade12age - age of youngest clade of size >= 12                ####
####                 (proxy for partial sweeps)                         ####
####  - <pop>_clade12age - pop-specific clade12age                      ####
####                       (proxy for pop-secific sweep                 ####
####  - <pop>_enrich - highest enrichment score for pop of any clade    ####
####  -                ( - log10(pval) under hypergeometric null)       ####      
#### Notes:                                                             ####
####  * all clade12age values are normalized by clade60age              ####
####  * hypergeometric p-values are for p(count >= observed count)      ####
############################################################################
library("ape", "phytools")
library("castor")
############################################################################
# get dir of this script (and also the treeStatFunctions.R script)
# - different code for interactive or batch (Rscript)
script_dir <- "./scripts/"
if(interactive()) {
  script_dir <- dirname(parent.frame(2)$ofile)
} else {
  library(tidyverse)
  get_script_dir <- function(){
    commandArgs() %>% 
       tibble::enframe(name=NULL) %>%
       tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
       dplyr::filter(key == "--file") %>%
       dplyr::pull(value) %>%
       dirname()
  }
  script_dir <- get_script_dir()
}
# print(script_dir)
############################################################################

############################################################################
# script arguments and i/o
# load functions
source(paste(script_dir,"treeStatFunctions.R",sep="/"))
args         <- commandArgs(trailingOnly = TRUE) # TRUE
# example: args         <- c("argTreeFiles/Contig1.50001-1650000.1000.tre.gz","argStats")
treeFile     <- args[1]
statsDir     <- args[2]
############################################################################

############################################################################
# set these based on your specific data set
total_n    <- 120      # total number of haploid samples in data set
species_n  <- 24       # number of haploid samples per species
############################################################################


############################################################################
# output stats file
statsFile <- sub(".*/Contig", "Contig", treeFile)
statsFile <- sub("\\.tre(\\.gz)?$", ".stats", statsFile)
statsFile <- paste(statsDir,statsFile,sep="/")
dir.create(statsDir,showWarnings=FALSE)
res <- file.create(statsFile,force=TRUE)
############################################################################

############################################################################
# other setup
# pop-individual key: note that pop key is already modified and pruned
pop_ind_key  <- read.csv("infoTables/individual-species-key-modified.txt",sep="\t")
popNames     <- unique(as.vector(pop_ind_key[,"Species"]))
scaffold     <- sub(".*/Contig", "Contig", treeFile)
scaffold     <- sub("\\..*", "", scaffold)
trees        <- ape::read.tree(treeFile)
coordinates  <- names(trees)
# write header
# OLD NAMES: treeStats    <- c("chrom", "pos", "clade60age", "clade12age_norm", paste(popNames,"clade12age_norm",sep="_"), paste(popNames, "enrich", sep="_"))
treeStats    <- c("chrom", "pos", "TMRCAH_all", "RT12", paste(popNames,"RTH",sep="_"), paste(popNames, "enrich", sep="_"))
write.table(t(treeStats), file=statsFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
############################################################################


############################################################################
# main loop - read all trees
counter=1
for(coord in coordinates) {
  # replace individual labels with pop labels and prune individuals not in key
  tree <- trees[[coord]]
  tree <- switchToPopLabels(tree , pop_ind_key)
  totalTreeCounts <- getSubtreeStats(tree, popNames)

  tree_stats_df <- data.frame(matrix(nrow=1,ncol=length(treeStats)))
  colnames(tree_stats_df) <- treeStats
  tree_stats_df$chrom  <- scaffold
  tree_stats_df$pos    <- coord

  # secondary loop - get subtree stats
  subtrees <- ape::subtrees(tree)
  subtree_stats_df <- getSubtreeStats(NULL,popNames)
  for(subtree in subtrees) {
    subtree_stats_df <- rbind( subtree_stats_df , getSubtreeStats(subtree,popNames) )
  }
  # replace NAs with 0s
  subtree_stats_df[is.na(subtree_stats_df)] <- 0
  # add enrichment scores
  subtree_stats_df <- computePopEnrichment(subtree_stats_df, popNames, totalTreeCounts) 

  # TMRCAH_all and RT12
  tree_stats_df$TMRCAH_all <- getCladeNage(subtree_stats_df,"total",n=total_n/2)
  tree_stats_df$RT12 <- round(getCladeNage(subtree_stats_df,"total",n=species_n/2) / tree_stats_df$TMRCAH_all , digits=4)
  # RTH' and enrichment scores
  for(pop in popNames) {
    pop_RTH <- paste(pop,"RTH",sep="_")
    pop_enrich     <- paste(pop,"enrich",sep="_")
    tree_stats_df[1,pop_RTH] <- round(getCladeNage(subtree_stats_df,pop,n=species_n/2) / tree_stats_df$TMRCAH_all , digits=4)
    tree_stats_df[1,pop_enrich]     <- max(subtree_stats_df[,pop_enrich])
  }

  # write stats line
  write.table(tree_stats_df[1,treeStats], file=statsFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  counter = counter+1
} # end of for(coord)
############################################################################

# gzip file
system(paste("gzip",statsFile,"-f"))

############################################################################

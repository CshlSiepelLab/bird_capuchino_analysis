### R script that contains various functions for computing cross coalescence statistics  ###
###  - the functions ape::mrca and getTreeDepth are used in ages calculations            ###  
library("ape", "phytools")

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
# load functions
source(paste(script_dir,"../smc_to_stats/treeStatFunctions.R",sep="/"))


###########################################################################################
#### getRecentCrossCoals - computes a given number of youngest cross coalescence events
####     between two populations.
####     Returns the num_events youngest coalescence ages between two populations
getRecentCrossCoals <- function(tree, pop1, pop2, num_events){
  # get MRCAs for every pair of tips in the tree (ape::mrca)  
  mrcas        <- ape::mrca(tree, full = FALSE)
  # subset rows and columns based on individuals from the two populations
  mrcas_subset <- mrcas[  which(rownames(mrcas)==pop1) ,  which(colnames(mrcas)==pop2) ] 
  # get table of unique MRCA ids (col mrcas_subset) and their counts (col Freq)
  mrcas_list   <- as.data.frame(table(mrcas_subset))
    
  # calculate ages of the MRCAs
  for(l in 1:nrow(mrcas_list)){
    mrcas_list[l,"age"] <- getTreeDepth(tree, idx=as.numeric(as.character(mrcas_list[l,1])), is.node=TRUE)
  }

  # create version of MRCA ages with repetitions
  mrca_ages <- rep(mrcas_list$age, mrcas_list$Freq)

  # returned first num_events after sorting from young to old
  return (sort(mrca_ages)[1:num_events])
}
###########################################################################################

###########################################################################################
#### getCrossCoalDistribution - calculates CC distribution per list of trees,
####   - Calls getRecentCrossCoals on each tree
####   - normalizes ages by clade60age, if required
####   - aggregates across trees in list
getCrossCoalDistribution <- function (trees, pop1, pop2, pop_ind_key, num_events, normalize = TRUE){
   cross_coals <- c()
   popNames    <- unique(as.vector(pop_ind_key[,"Species"]))
   for (t in 1:length(trees) ){
      tree <- switchToPopLabels(trees[[t]], pop_ind_key )
      youngest_cross_coals <- getRecentCrossCoals(tree, pop1, pop2, num_events )
      if(normalize == TRUE){
         # calaulate TMRCAH of entire tree by getting clades and finding
         # youngest one that contains at least half the leaves
         subtrees <- ape::subtrees(tree)
         subtree_stats_df <- getSubtreeStats(NULL,popNames)
         for(subtree in subtrees) {
            subtree_stats_df <- rbind( subtree_stats_df , getSubtreeStats(subtree,popNames) )
         }
         # replace NAs with 0s
         subtree_stats_df[is.na(subtree_stats_df)] <- 0
         tmrcah                <- getCladeNage(subtree_stats_df,"total",n=length(tree$tip.label)/2)
         youngest_cross_coals  <- youngest_cross_coals / tmrcah
      }
      # aggregate cross coalescence ages
      cross_coals <- c(cross_coals, youngest_cross_coals)
  }

  return(cross_coals)
}
###########################################################################################

###########################################################################################
#### computeQuantileDiff - computes quantile difference between two vectors of values
####    quantile difference is the absolute value of the difference between 1/2 and
####    the quantile that the median value of one distribution achieves in the second.
####    This difference is computed in both directions and the maximum value is returned
####    Returned value is set such that sign is + if first list is greater than second
computeQuantileDiff <- function(valueList1, valueList2){

   # get median and quantile
   length1        <- length(valueList1)
   length2        <- length(valueList2)
   median1        <- median(valueList1)
   median2        <- median(valueList2)
   quantile1      <- length(which(valueList1 >= median2))/length1
   quantile2      <- length(which(valueList2 < median1))/length2
   if(min(quantile1 , quantile2)<0.5) {
      return(min(quantile1 , quantile2) - 0.5)
   }
   return (max(quantile1 , quantile2) - 0.5)
}
###########################################################################################



###########################################################################################
#### getTreesPerRegion - returns list of trees per given genomic region
####   - Aggregates relevant trees across different tree files
####   - normalizes ages by clade60age, if required
####   - aggregates across trees in list
getTreesPerRegion <- function (treeDir, scaffold, startPos, endPos) {
   # get all tree files for scaffold
   # tree files names: <scaffoldID>.<start>-<end>.*
   # extract also coordinates
   treeFiles       <- list.files(treeDir, pattern = paste(sep='',"^",scaffold,"[.]" ))
   treeFilesStart  <- as.numeric(sub(".*[.]", "", (sub("-.*", "", treeFiles))) ) - 1
   fileOrder       <- order(treeFilesStart)
   treeFiles       <- treeFiles[fileOrder]
   treeFilesStart  <- treeFilesStart[fileOrder]
   treeFilesEnd    <- as.numeric(sub("[.].*", "", (sub(".*-", "", treeFiles))) )
   relevantInds    <- which(treeFilesStart <= endPos & treeFilesEnd >= startPos)
   treeFiles       <- treeFiles[relevantInds]
   trees <- c()
   for (file in treeFiles){
      trees_from_file <- ape::read.tree(paste(sep='',treeDir,"/",file))
      trees_from_file <- trees_from_file[startPos <= as.numeric(names(trees_from_file)) & endPos > as.numeric(names(trees_from_file))]
      trees <- c(trees , trees_from_file)
   }
   return(trees)
}
###########################################################################################



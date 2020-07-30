### R script that contains various functions for computing tree and subtree statistics  ###
library("ape", "phytools")
library("castor")
library("plyr")



##################################################################################
####  switchToPopLabels  - returns tree with tip labels replaced by population labels
####                       also prunes all tips that do not belong to populations in key
switchToPopLabels <- function(tree, pop_ind_key) {
  # replace tip labels with population labels for both haplotypes of every individual
  popNames <- unique(pop_ind_key$Species)
  tree$tip.label <- mapvalues(tree$tip.label,as.vector(pop_ind_key$Individual_hap),as.vector(pop_ind_key$Species),warn_missing=FALSE)
  # prune tips that do not belong to pops
  non_pop_tips <- which(is.na(match(tree$tip.label,popNames)))
  return(drop.tip(tree,non_pop_tips))
}
##################################################################################

##################################################################################
####  getTreeDepth  - return TMRCA for given tree - function from argweaver  ####
getTreeDepth <- function(tree, idx=-1, is.node=FALSE) {
  if (class(tree) != "phylo") {
    if (!requireNamespace("ape", quietly=TRUE)) {
      stop("ape package is required for getTreeDepth")
    }
    tree <- ape::read.tree(text=tree)
  }
  if (idx==-1) {
    # find root
    node <- unique(tree$edge[,1][which(!is.element(tree$edge[,1], tree$edge[,2]))])
  } else if(is.node == FALSE){
    # find given node
    node <- tree$edge[idx,2]
  } else {
    node <- idx
  }
  anc <- which(tree$edge[,1]==node)
  # if not leaf
  if (length(anc)==2) {
    return(max(tree$edge.length[anc[1]] + getTreeDepth(tree, anc[1]),
               tree$edge.length[anc[2]] + getTreeDepth(tree, anc[2])))
  }
  0
}
##################################################################################

##################################################################################
####  getSubtreeStats - computes stats for subtree: leaf counts per pop, total leaf count, and age
getSubtreeStats <- function (tree,popNames){
  popCounts <- paste(popNames,"count",sep="_")
  stats <- c(popCounts, "total_count", "age");
  if(is.null(tree)) {
    subtree_stats_df <- data.frame(matrix(ncol= length(stats), nrow=0))
    colnames(subtree_stats_df) <- stats
  } else {
    subtree_stats_df <- data.frame(matrix(ncol= length(stats), nrow=1))
    colnames(subtree_stats_df) <- stats
    pop_counts <- table(tree$tip.label)
    subtree_stats_df[1,popCounts]     <- pop_counts[popNames]
    subtree_stats_df[1,"total_count"] <- length(tree$tip.label) 
    subtree_stats_df[1,"age"]         <- getTreeDepth(tree, idx=-1 ) 
  }
  return(subtree_stats_df);
}
##################################################################################

##################################################################################
####  getCladeNage - returns age of youngest clade that contains n individuals of certain type
##     receives a data frame with subtree ages and counts and a column with relevant counts
getCladeNage <- function (subtree_stats_df,countType,n){
  countColumn <- paste(countType,"count",sep="_")
  cladeNsubtrees <- which(subtree_stats_df[,countColumn] >= n)
  return (min(subtree_stats_df[cladeNsubtrees,"age"]))
}
##################################################################################

##################################################################################
####  computePopEnrichment - computes population enrichment score
##     receives a data frame with subtree counts for all population, and relevant counts
##     enrichment score is minus log10(p-val) under hypergeometric distribution
computePopEnrichment <- function (subtree_stats_df,popNames, sampleCounts){
  # compute enrichment score for each pop as minus log10(p-val) under hypergenometric null distribution
  subtree_total_counts <- subtree_stats_df[,"total_count"]
  
  for(pop in popNames) {
    count_col  <- paste(pop,"count",sep="_")
    subtree_pop_counts   <- subtree_stats_df[,count_col]

    # compute P(count >= observed count)
    subtree_stats_df[,paste(pop,"enrich",sep="_")] <- 
       dhyper(x=subtree_pop_counts, 
              m=sampleCounts[1,count_col],
              n=sampleCounts[1,"total_count"]-sampleCounts[1,count_col],
              k=subtree_total_counts) +
       phyper(q=subtree_pop_counts, 
              m=sampleCounts[1,count_col],
              n=sampleCounts[1,"total_count"]-sampleCounts[1,count_col],
              k=subtree_total_counts,
              lower.tail=FALSE)
    subtree_stats_df[,paste(pop,"enrich",sep="_")] <- round(digits=3, (-log10(subtree_stats_df[,paste(pop,"enrich",sep="_")]))) 

  }

  return(subtree_stats_df)

}
##################################################################################


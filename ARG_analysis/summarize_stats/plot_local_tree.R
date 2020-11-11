#############################################################################################
#### - Extract and plot significant trees for a given range                              ####
#### - Every 20 kb block in the given range is analyzed separately                       ####
#### - In case that more than one population has significant values of enrichment        ####
####      and/or RTH, several trees will be printed for that 20kb region.                ####
#### - The output - a single pdf file with the most significant trees found in the range ####
#############################################################################################
library("ape", "phytools", "argweaver")
# also requires RPHAST for subsetting tree.
# library("gridExtra", "grid")


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
# load functions for getting trees
source(paste(script_dir,"../smc_to_stats/crossCoalFunctions.R",sep="/"))


#############################################################################################
####                                    SCRIPT ARGUMENTS                                 ####
#############################################################################################
args <- commandArgs(trailingOnly = TRUE)
# args <- c("Contig252" , "420000" , "510000", "plots/treeExample.pdf") #EXAMPLE ARGUMENT LINE FOR PEAK ON Contig252
scaffold     <- args[1]
startPos     <- as.numeric(args[2])
endPos       <- as.numeric(args[3])  
outFile      <- args[4]

#############################################################################################
####                                         SETUP                                       ####
#############################################################################################

statFileDir   <- "./argStats/"
treeFileDir   <- "./argTreeFiles/"

# compute list to map individuals to colors.
pop_ind_key   <- read.csv("infoTables/individual-species-key-modified.txt",sep="\t", stringsAsFactors = FALSE)
populations   <- unique(pop_ind_key[,"Species"])
ind_color_list <- list()
for (i in 1:nrow(pop_ind_key)) {
   ind_color_list[[pop_ind_key$Individual_hap[i]]] <- pop_ind_key$Color[i]
}

# Change here if you wish to prune tree for a specific subset of individuals (NULL indicates all individuals)
inds_to_keep <- NULL

#############################################################################################
####                                 FUNCTIONS USED IN SCRIPT                            ####
#############################################################################################


###########################################################################################
#### myPlotTree - receives a tree in text from tree file and plots id
##      prints the tree via argweaver::plotTree
##      additional details are plotted at the bottom  
myPlotTree <- function(tree, ind_color_list, tips_to_keep=NULL,
                              additional_text="", cex.leafname=0.3, lwd=2, mar=c(4,4,1,1)){
                              
   if(is.null(tips_to_keep)) {
      tips_to_keep <- names(ind_color_list)
   }

   argweaver::plotTree(tree=tree, col=ind_color_list, keepSeqs=tips_to_keep, lwd=lwd, cex.leafname=cex.leafname, mar=mar, leafLabels=NULL) 
    
   if(additional_text != "") {
      mtext(additional_text, side=1, line=2, adj=0.5,cex = 1.0)
   }

}


###########################################################################################
#### getTreesAndStatsPerRegion - returns list of trees per given genomic region 
####                                 + a list of corresponding stats
####   - Aggregates relevant trees across different tree files
####   - Does the same for stat files
####   - returns a list with two items: trees - list of trees ; stats - table of stats 
####   - Derived from function getTreesPerRegion from crossCoalfunctions in directory smc_to_stats
getTreesAndStatsPerRegion <- function (treeDir, statDir, scaffold, startPos, endPos) {
   # get all tree files for scaffold
   # tree files names: <scaffoldID>.<start>-<end>.*
   # extract also coordinates
   treeFiles       <- list.files(treeDir, pattern = paste(sep='',"^",scaffold,"[.]" ))
   statFiles       <- list.files(statDir, pattern = paste(sep='',"^",scaffold,"[.]" ))
   treeFilesStart  <- as.numeric(sub(".*[.]", "", (sub("-.*", "", treeFiles))) ) - 1
   fileOrder       <- order(treeFilesStart)
   treeFiles       <- treeFiles[fileOrder]
   treeFilesStart  <- treeFilesStart[fileOrder]
   treeFilesEnd    <- as.numeric(sub("[.].*", "", (sub(".*-", "", treeFiles))) )
   relevantInds    <- which(treeFilesStart <= endPos & treeFilesEnd >= startPos)
   treeFiles       <- treeFiles[relevantInds]
   statFiles       <- statFiles[relevantInds]
   # print(statFiles)
   trees <- c()
   for (file in treeFiles){
      # DO NOT USE ape::read.tree TO RETAIN ARG-BASED INFORMATION FROM TREE
      trees_from_file        <- read.table(paste(sep='',treeDir,"/",file),header=FALSE, stringsAsFactors=FALSE, quote="")
      names(trees_from_file) <- c("pos","tree")
      trees_from_file <- trees_from_file[startPos <= trees_from_file$pos & endPos > trees_from_file$pos,]
      trees <- c(trees , trees_from_file)
   }
   stats <- c()
   for (file in statFiles){
      stats_from_file <- read.table(paste(sep='',statDir,"/",file),header=TRUE, stringsAsFactors=FALSE, quote="")
      stats_from_file <- stats_from_file[which(startPos <= stats_from_file$pos & endPos > stats_from_file$pos),]
      stats <- rbind(stats, stats_from_file)
   }
   result <- list("trees" = trees, "stats" = stats)
   return(result)
}
###########################################################################################


###########################################################################################
#### getBestTrees - returns a data frame with trees best demonstrating
####                 high species enrichment and low RTH
getBestTrees <- function (trees, stats, populations) {
   attributes<-c("pos","stat_val","tree")
   enrich_stats  <- paste(sep="_",populations,"enrich")
   rth_stats     <- paste(sep="_",populations,"RTH")
   stat_names    <- c(enrich_stats , rth_stats)
   trees_df<-data.frame(matrix(nrow=length(stat_names),ncol=length(attributes)))
   colnames(trees_df)<-attributes
   rownames(trees_df)<-stat_names
   
   # find "best" tree for each stat
   for(stat in stat_names) {
      # find MAX for enrichment scores / minimum for RTH
      isEnrich <- length(grep("_enrich$",stat)) > 0
      stat_row <- order(stats[,stat],decreasing=isEnrich)[1]
      pos      <- stats$pos[stat_row]
      trees_df[stat,"pos"]        <-  pos
      trees_df[stat,"stat_val"]   <-  stats[stat_row,stat]
      ## cat("best tree for stat ",stat," found in row ",stat_row," position ",pos," value ", stats[stat_row,stat],"\n")
      tree_row <- which(trees$pos==pos)
      if(length(tree_row) != 1) {
         trees_df[stat,"tree"]   <-  NULL
      } else {
         ## cat(paste(stat, stats[stat_row,stat], pos,stat_row,tree_row,"\n"))
         trees_df[stat,"tree"]   <-  trees$tree[tree_row]
      }      
   }
   # order trees based on decreasing enrichment scores.
   # for each population show most enriched tree and then lowest RTH tree
   row_order <- 1:length(stat_names)
   row_order[2*1:length(populations)-1]  <- order(trees_df[enrich_stats,"stat_val"],decreasing=TRUE)
   row_order[2*1:length(populations)]  <- row_order[2*1:length(populations)-1] + length(populations)
   trees_df  <-  trees_df[row_order,]
      

  return(trees_df)
}


#############################################################################################
####                                     MAIN BODY                                       ####
#############################################################################################

cat("Loading trees and stats for region: ",scaffold,":",startPos,"-",endPos,"\n")
results  <- getTreesAndStatsPerRegion (treeFileDir, statFileDir, scaffold, startPos, endPos)
cat("Extracting best demonstrative tree for each stat\n")
trees_df <- getBestTrees (results$trees , results$stats, populations)

cat("Generating plot for every tree:\n")
pdf(outFile, height=5, width=nrow(pop_ind_key)/10)
for(i in 1:nrow(trees_df)) {
   caption <- paste(sep="",scaffold,":",trees_df$pos[i],"   (",rownames(trees_df)[i]," = ", trees_df$stat_val[i],")")
   cat(" - ", caption,"\n")
   myPlotTree(trees_df$tree[i], ind_color_list, additional_text=caption, tips_to_keep=inds_to_keep)
}
dev.off()


##########################################################################################################
###################  summarize cross coalescnece differences in Fst peaks and control blocks   ###########
##########################################################################################################

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
source(paste(script_dir,"../smc_to_stats/crossCoalFunctions.R",sep="/"))

args     <- commandArgs(trailingOnly = TRUE) # TRUE
treeDir          <- args[1]
regionID         <- args[2]
scaffold         <- args[3]
leftFlankEnd     <- as.numeric(args[4])
rightFlankStart  <- as.numeric(args[5])
startPos         <- as.numeric(args[6])
endPos           <- as.numeric(args[7])
outFile          <- args[8]
flankLen                       <- 200000

pop_ind_key  <- read.csv("infoTables/individual-species-key-modified.txt",sep="\t")
popNames     <- unique(as.vector(pop_ind_key[,"Species"]))
popPairs     <- matrix("",nrow=length(popNames)*(length(popNames)-1)/2,ncol=2)
k            <- 1
for(i in 1:(length(popNames)-1)) {
   for(j in (i+1):length(popNames)) {
      popPairs[k,] <- popNames[c(i,j)]
      k <- k+1
   }
}

df_cols           <- paste(sep="_",popPairs[,1],popPairs[,2])
out_df            <- data.frame(matrix(nrow=1,ncol=length(df_cols)))
colnames(out_df)  <- df_cols
rownames(out_df)  <- regionID
cat("summarizing cross coals for region ",regionID,": ")
cat(startPos,"-",endPos)

trees            <- getTreesPerRegion(treeDir, scaffold, leftFlankEnd-flankLen, rightFlankStart+flankLen)
coords           <- as.numeric(names(trees))
minCoord         <- min(coords)-(coords[2]-coords[1])/2
maxCoord         <- max(coords)+(coords[2]-coords[1])/2
leftFlankEnd     <- max(minCoord, leftFlankEnd)
leftFlankStart   <- max(minCoord, leftFlankEnd-flankLen/2)
rightFlankStart  <- min(maxCoord, rightFlankStart)
rightFlankEnd    <- min(maxCoord, rightFlankStart+flankLen/2)
# cat("Flanking region for ",peakID,": ",leftFlankStart,"-",leftFlankEnd,"   and   ",rightFlankStart,"-",rightFlankEnd,"\n")
if(leftFlankEnd-leftFlankStart < flankLen/2) {
   rightFlankEnd <- min(maxCoord,rightFlankStart+flankLen-(leftFlankEnd-leftFlankStart))
   # cat("elong right ",rightFlankEnd,"\n")
}
if(rightFlankEnd-rightFlankStart < flankLen/2) {
   leftFlankStart <- max(minCoord,leftFlankEnd-flankLen+(rightFlankEnd-rightFlankStart))
   # cat("elong left ",leftFlankstart,"\n")
}
if(leftFlankEnd-leftFlankStart + rightFlankEnd-rightFlankStart < flankLen) {
   cat(" Flanking regions ",leftFlankStart,"-",leftFlankEnd,"   and   ",rightFlankStart,"-",rightFlankEnd," not long enough\n")
} else {
   cat(" Flanking regions ",leftFlankStart,"-",leftFlankEnd,"   and   ",rightFlankStart,"-",rightFlankEnd)
   peakTrees    <- trees[which(coords <= endPos & coords > startPos)]
   flankTrees   <- trees[c( which(coords <= leftFlankEnd & coords > leftFlankStart) , which(coords <= rightFlankEnd & coords > rightFlankStart) )]
   cat(" ",length(peakTrees)," peak trees and ",length(flankTrees)," flank trees")
   for(p in 1:nrow(popPairs)) {
      peakCCs   <- getCrossCoalDistribution(peakTrees , popPairs[p,1], popPairs[p,2], pop_ind_key, num_events=10, normalize = TRUE)
      flankCCs  <- getCrossCoalDistribution(flankTrees, popPairs[p,1], popPairs[p,2], pop_ind_key, num_events=10, normalize = TRUE)
      perDiff   <- computeQuantileDiff(peakCCs , flankCCs)
      out_df[1,paste(sep="_",popPairs[p,1], popPairs[p,2])] <- round(perDiff,digits=2)
      cat("\n  ",paste(sep="_",popPairs[p,1], popPairs[p,2])," : ",perDiff)
   } # end of for(p)
   cat("\n")
} # end of if(FLANKING REGION NOT LONG ENOUGH)
write.table(out_df, outFile, quote=FALSE, sep="\t")

warnings()


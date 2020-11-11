  ### R script to run on Linux  ###
############################################################################
##### compute average cross coal stats per window across given region   ####
#####  - reads off of tree files                                        ####
#####  - averages stats across windows of given length (20 kb)          ####
#####  - computes start and end index of windows based on start and end ####
#####    indices of tree file                                           ####
#####  - produces BED file with average stats per window                ####
#####  - checks that all windows are covered equally by sampled trees   ####
#####  - writes information on scaffold coverage in end of given info file #
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
source(paste(script_dir,"crossCoalFunctions.R",sep="/"))


############################################################################
# script arguments and i/o
# example:  args         <- c("./argTreeFiles","Contig252", "220000", "710000", "20000", "argStats-windowed")
args <- commandArgs(trailingOnly = TRUE) # TRUE
treeDir     <- args[1]
scaffold    <- args[2]
startPos    <- as.numeric(args[3])
endPos      <- as.numeric(args[4])
window_len  <- as.numeric(args[5])
outDir      <- args[6]
# rounding parameter
roundDigits  <- 3   # keep 3 decimal digits in age-based score
############################################################################

############################################################################
# determine all population pairs
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

############################################################################
# determine start index of windows based on start index of first tree file in
# scaffold plus integer multiples of window length
treeFiles       <- list.files(treeDir, pattern = paste(sep='',"^",scaffold,"[.]" ))
treeFilesStart  <- as.numeric(sub(".*[.]", "", (sub("-.*", "", treeFiles))) ) - 1
scaffoldStart   <- min(treeFilesStart)
windowStart     <- scaffoldStart + window_len * floor((startPos - scaffoldStart)/window_len)
windowEnd       <- windowStart   + window_len * ceiling((endPos - windowStart)/window_len)

############################################################################
# create output file
outFile         <- paste(outDir,"/",scaffold,".",windowStart,"-",windowEnd,".crosscoal.bed",sep="")
stat_list       <- paste(sep="_",popPairs[,1],popPairs[,2])
window_stats    <- c("chrom", "startPos", "endPos", stat_list)
dir.create(outDir,showWarnings=FALSE)
file.create(outFile,showWarnings=FALSE)
write.table(t(window_stats), file=outFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)

############################################################################
# get list of relevant trees
trees             <- getTreesPerRegion(treeDir, scaffold, windowStart, windowEnd)
coords            <- as.numeric(names(trees))
# expected number of trees per window based on gap between first two trees
num_trees_per_seg <- window_len / (coords[2] - coords[1])



############################################################################
# main loop - iterate over windows
# pos holds current starting position for window
pos <- windowStart
while(pos < windowEnd) {
   window_lines <- which(coords > pos & coords <= pos+window_len)
   if(length(window_lines) == 0) {
      pos = pos + window_len
      next
   } else if(length(window_lines) != num_trees_per_seg) {
      # unexpected number of stat lines in window
      message <- paste(scaffold,"unexpected number of trees (",length(window_lines)," instead of",num_trees_per_seg,") in window strating at",pos,".\n")
      cat(message)
      break
   }
   # fill window stats - initialize chrom and positions
   window_stats_df <- data.frame(matrix(nrow=1,ncol=length(window_stats)))
   colnames(window_stats_df)    <- window_stats
   window_stats_df$chrom        <- scaffold
   window_stats_df$startPos     <- pos
   window_stats_df$endPos       <- pos+window_len
   # compute mean cross coal across trees in window per pop pair
   for(p in 1:nrow(popPairs)) {
      # cat("Calculating CCs in window ",scaffold,":",pos,"-",pos+window_len,"for pair ",popPairs[p,1],"-",popPairs[p,2],"\n")
      CCs   <- getCrossCoalDistribution(trees[window_lines] , popPairs[p,1], popPairs[p,2], pop_ind_key, num_events=10, normalize = TRUE)
      window_stats_df[1,paste(sep="_",popPairs[p,1], popPairs[p,2])] <- round(mean(CCs),roundDigits)
      # print(window_stats_df)
   } # end of for(p)
    
   # write stats line and advance pos
   write.table(window_stats_df[1,window_stats], file=outFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
   pos = pos + window_len
} # end of while(pos)

if(pos > windowEnd) {
   # unexpected gap in the end
   message <- paste(scaffold," gap in end of region",i,"( pos",windowEnd,"-",pos,")\n")
   cat(message)
}
############################################################################

# gzip file - TURNED OFF
# system(paste("gzip",outFile,"-f"))
############################################################################

## computes a stats file for every trees file
# writes commands in bash file and groups according to batch size, ordered by block length to balance load across parallel runs
#
# command line example: Rscript scripts/generateTreeStats.R argTreeFiles/Contig0.50001-1850000.1000.tre.gz argStats

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

outFile     <- "./create_stat_files.sh"
script_path <- paste(script_dir,"tre_to_stats.R",sep="/");
tree_dir    <- "./argTreeFiles"
stats_dir    <- "./argStats"
treeFiles   <- list.files(tree_dir,"*.tre(.gz)?$")
# compute table with tree file and segment length 
coords      <- sub("[^.]*\\.","",treeFiles)
coords      <- sub("\\..*","", coords)
startInd    <- as.numeric(sub("-[0-9]*","",coords))
endInd      <- as.numeric(sub("[0-9]*-","",coords))
table       <- cbind(treeFiles, endInd-startInd+1)
colnames(table) <- c("treeFile", "segmentLength")
table       <- table[order(as.numeric(table[,"segmentLength"])),] 

# script parameters:
numJobs     <- 20         # num jobs to run in parallel

write("#!/bin/bash", file = outFile, append=FALSE)

for (i in  1:nrow(table)){ 
  command <-sprintf("Rscript %s %s/%s %s &", script_path, tree_dir, table[i,"treeFile"], stats_dir)
  write(command, file=outFile , append = TRUE)
  if (i %% numJobs == 0){
    # the "wait" means it waits until all previous commands are finished and then moves to the next command.
    write("wait", file=outFile, append = TRUE)
  }
}

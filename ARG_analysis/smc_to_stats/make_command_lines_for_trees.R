## read trimmed arg block list (ARGblock-coordinates-trimmed.txt) and create tree filtered tree files for each bed file
# writes commands in bash file and groups according to batch size, ordered by block length to balance load across parallel runs

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

inputFile   <- "./infoTables/ARGblock-coordinates-trimmed.txt"
outFile     <- "./create_tree_files.sh"
lines       <- read.table(inputFile)
script_path <- paste(script_dir,"bed_to_tre.R",sep="/");
bed_dir     <- "./argBedFiles"
tree_dir    <- "./argTreeFiles"

# script parameters: <-- SET THESE, IF NEEDED
skipInd     <- 500        # interval between sampled trees
numJobs     <- 20         # num jobs to run in parallel

# the input file:
lines <- read.table(inputFile)
tableColumns<- c("file", "scaffold", "startInd", "endInd", "blockLen"); 
lines$V5 <- as.numeric(lines$V4) - as.numeric(lines$V3) + 1
colnames(lines)<-tableColumns
lines$file <- sub(".*/",paste(bed_dir,"/",sep=""),lines$file);
lines$file <- sub(".smc.gz$",".bed.gz",lines$file);

# sort by block length for better balancing between parallel jobs, and filter blocks with size 0 (all ARG block trimmed)
lines <- lines[order(lines$blockLen),]
lines <- lines[which(lines$blockLen>0),]
write("#!/bin/bash", file = outFile, append=FALSE)

for (i in  1:nrow(lines)){ 
  command <-sprintf("Rscript %s %s %s %s %s %s &", script_path, lines[i,"file"], tree_dir, lines[i,"startInd"], lines[i,"endInd"], skipInd)
  write(command, file=outFile , append = TRUE)
  if (i %% numJobs == 0){
    # the "wait" means it waits until all previous commands are finished and then moves to the next command.
    write("wait", file=outFile, append = TRUE)
  }
}

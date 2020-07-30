## read args file list. (e.g. ARGblock-coordinates.txt or ARGblock-coordinates-trimmed.txt)
# For every row, get the contig name and iteration number and write a command line to a bash file to pass later to "call_smc_to_bed.R"
# potentially filter based on MCMC iteration and arg block length
# list all filtered arg blocks in separate file
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



inputFile   <- "./infoTables/ARGblock-coordinates.txt"
outFile     <- "./create_bed_files.sh"
infoFile    <- "./infoTables/arg_blocks_with_no_bed.txt"
lines       <- read.table(inputFile)
script_path <- paste(script_dir,"smc_to_bed.R",sep="/");
smc_dir     <- # <-- SET SOURCE DIR CONTAINING SMC FILES
log_dir     <- # <-- SET SOURCE DIR CONTAINING LOG FILES
out_dir     <- "./argBedFiles"

# script parameters:
numJobs        <- 20         # num jobs to run in parallel
mcmc_iter      <- 1000       # MCMC iteration for filtering
min_block_len  <- 100000     # minimum scaffold length for filtering


# the input file:
lines <- read.table(inputFile)
tableColumns<- c("file", "scaffold", "startInd", "endInd", "mcmcIter", "blockName", "blockLen"); 
lines$V5 <- sub(".*out.","",lines$V1)
lines$V5 <- as.numeric(sub(".smc.gz","",lines$V5))
lines$V6 <- sub(".*/","",lines$V1)
lines$V6 <- sub("_.*","",lines$V6)
lines$V7 <- as.numeric(lines$V4) - as.numeric(lines$V3) + 1
colnames(lines)<-tableColumns
# sort by block length for better balancing between parallel jobs
lines <- lines[order(lines$blockLen),]
write("#!/bin/bash", file = outFile, append=FALSE)
write(file=infoFile, append=FALSE, paste(sep="\t","file","mcmcIter","blockLen"))

for (i in  1:nrow(lines)){ 
  if(lines[i,"mcmcIter"] < mcmc_iter || lines[i,"blockLen"] < min_block_len) {
     write(file=infoFile, append=TRUE, paste(sep="\t",lines[i,"file"], lines[i,"mcmcIter"], lines[i,"blockLen"]))
  } else {
    command <-sprintf("Rscript %s %s %s %s 0 %s %s %s &", script_path, lines[i,"blockName"], lines[i,"mcmcIter"], lines[i,"mcmcIter"], smc_dir, log_dir, out_dir)
    write(command, file=outFile , append = TRUE)
    if (i %% numJobs == 0){
      # the "wait" means it waits until all previous commands are finished and then moves to the next command.
      write("wait", file=outFile, append = TRUE)
    }
  }
}

######################################################################################
## computes a window-stats file for all scaffolds
#  loops over all scaffolds and calls window_stats.R
######################################################################################

# get dir of this script
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


######################################################################################
# arguments
statsDir    <- "./argStats/"
outDir      <- "./argStats-windowed/"  
infoFile    <- "./infoTables/scaffold-window-coverage.txt"

# script parameters:
window_len <- 20000     # length of non-overlapping windows for averaging

# compute list of scaffolds and sort by number
statFiles   <- list.files(statsDir,"*.stats(.gz)?$")
scaffolds   <- sub(".*/","", statFiles)
scaffolds   <- sub("\\..*","", statFiles)
scaffolds   <- unique(scaffolds)
scaffoldIDs <- as.numeric(sub("[^0-9]*","",scaffolds))
scaffolds   <- scaffolds[order(scaffoldIDs)]

# main loop by scaffold
write("# information on scaffold blocks (missing and inconsistent blocking)", infoFile, append=FALSE)
for(scaffold in scaffolds) {
  args <- c(statsDir, scaffold, window_len, outDir, infoFile)
  source(paste(script_dir,"window_stats.R", sep="/"))
}

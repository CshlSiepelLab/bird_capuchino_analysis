### R script to run on Linux  ###
################################################################
#### based on make_tre_file_by_coordinates.R                ####
#### converts an ARG bed file into a tre file               ####
#### restructs to specific region (after trimming)          ####
#### subsamples tree after given skip                       ####
################################################################


################################################################
## inpit arguments
args <- commandArgs(trailingOnly = TRUE) # TRUE

bedFile   = args[1] #input bed file
outDir    = args[2] # out dir where all tree files are put
startInd  = args[3] # start index
endInd    = args[4] # end index
skipInd   = args[5] # skips between sampled trees
################################################################


################################################################
#      create tree file name from bed file name + inds
dir.create(outDir,showWarnings=FALSE)
mcmcIter <- sub(".*out.","",bedFile)
mcmcIter <- sub(".bed.gz","",mcmcIter)
scaffold <- sub(".*/","",bedFile)
scaffold <- sub("[.]bed[.]gz","",scaffold)
scaffold <- sub("_out\\..*","",scaffold)
scaffold <- sub("\\..*","",scaffold)
treeFile <- paste(sep="",outDir,"/",scaffold,".",startInd,"-",endInd,".",mcmcIter,".tre")
res      <- file.create(treeFile)
startInd <- as.numeric(startInd)
endInd   <- as.numeric(endInd)
skipInd  <- as.numeric(skipInd) # skips between sampled trees
# sampled indices start at (startInd-1)+(skipInd/2) and then skip every skipInd bases
indices <- seq(startInd-1+(skipInd/2), endInd, by = skipInd)
################################################################

################################################################
# single call to tabix and then loop over table
################################################################
  command <- sprintf("tabix %s %s:%.0f-%.0f", bedFile, scaffold, startInd, endInd)
  # print(command)
  treeTable <- read.table(pipe(command), header=FALSE, stringsAsFactors=FALSE)
  names(treeTable) <- c("chrom", "chromStart", "chromEnd", "MCMC", "tree")
  ind = startInd-1 + skipInd/2
  for (i in 1:nrow(treeTable)){
    while(ind <= treeTable[i,"chromEnd"] && ind <= endInd ) {
       write(paste(format(ind,scientific=FALSE) , treeTable[i,"tree"],sep="\t"), file= treeFile, append = TRUE)
       ind = ind + skipInd
    }
  }
################################################################

################################################################
# gzip tree file
system(paste("gzip",treeFile,"-f"));
################################################################

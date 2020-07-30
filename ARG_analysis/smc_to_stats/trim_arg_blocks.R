#################################################################################
##  re-mapping args regions ##
## read a txt file with regions for every ARG and create a new version of the file, with regions updated so that edges are trimmed
## and the total region divides to 20kb size segments.
## for args that are part of a large block, divide the block to regions without repeats

## args with iteration less than "iteration_number" are discarded.

## returns also a txt with summary of the total length, num.of ARGs in size less than 100kb, num. of ARGs in size more than 100kb 
# and num. of ARGs with mcmc iteration less than 1000
#################################################################################

# input / output file names:
workDir   <- "./infoTables"
inputFile <- paste(sep="/", workDir,"ARGblock-coordinates.txt")
outFile   <- paste(sep="/", workDir,"ARGblock-coordinates-trimmed.txt")
infoFile  <- paste(sep="/", workDir,"ARGblock-coordinates-info.txt")
lines     <- read.table(inputFile)
# script parameters:
mcmc_iter      <- 1000
min_block_len  <- 100000
segment_len    <- 20000
trim_len       <- 50000


divideBlocksToRegionsBySegments <- function (args_df, trim_len, segment_len){
  ## input: df with regions data on all args in the scaffold
  ## output: df with the regions data on all args in the scaffold, after cutting the and dividing to new regions
  ## works on a scaffold that contains one or more ARG blocks
  ## 0. trim left side of first ARG block in scaffold by  setting start index to trim_len
  ## 1. set end index of ARG block by tiling ARG block with segments of length segment_len
  ##    from start index to closest index to end index minus trim_len
  ## 2. set start index of next ARG block to end inxed of previous block and go to (1)
  ## If an ARG block is missing, then we reset process

  args_df$blockIndex <- 0;
  if(nrow(args_df) > 1) {
    args_df$blockIndex <- as.integer(sub(".*\\.","",args_df$blockName))
  }
  ordered_args_df <- args_df[order(args_df$blockIndex),]

  for(i in 1:nrow(ordered_args_df)){
  
    # If ARG block is first in scaffold, or after a missing block, 
    # then set start index accroding to trim_len.
    # Otherwise, set start index according to end index of previous block

    if(i==1 || ordered_args_df[i,"blockIndex"] != ordered_args_df[i-1,"blockIndex"] + 1 ) {
      ordered_args_df[i,"startInd"] = ordered_args_df[i,"startInd"] + trim_len;
    } else {
      ordered_args_df[i,"startInd"] = ordered_args_df[i-1,"endInd"]+1;
    }

    # set end index of ARG block by tiling segments of length segment_len
    num_segments <- round((ordered_args_df[i,"endInd"] - trim_len - ordered_args_df[i,"startInd"]) / segment_len)
    ordered_args_df[i,"endInd"] = ordered_args_df[i,"startInd"] +  num_segments * segment_len - 1 
  }
  return(ordered_args_df)
}
#################################################################################

#################################################################################
total_length_contigs <- function (lines_from_regions_file){
  # returns the total length of the ARGs

  total_length <- sum(lines_from_regions_file$endInd - lines_from_regions_file$startInd)

  return(total_length);
}
#################################################################################


#################################################################################
#                            MAIN SCRIPT CODE                                   #
#################################################################################


#add columns with iteration_number and ARG block name
tableColumns<- c("file", "scaffold", "startInd", "endInd", "mcmcIter", "blockName"); 
lines$V5 <- sub(".*out.","",lines$V1)
lines$V5 <- sub(".smc.gz","",lines$V5)
lines$V6 <- sub(".*/","",lines$V1)
lines$V6 <- sub("_.*","",lines$V6)
colnames(lines)<-tableColumns
num_scaffolds <- length(unique(lines$scaffold))

# filter based on iteration
lines_mcmc_iter     <- subset(lines, as.numeric(lines$mcmcIter) >= mcmc_iter)
lines_mcmc_iter_not <- subset(lines, as.numeric(lines$mcmcIter) <  mcmc_iter)
# get ARG block names and make sure they are unique
arg_block_names <- lines_mcmc_iter$blockName;


# Not completely sure this is needed. We're now filtering this out before processing the files.
arg_blocks_few_iterations <- lines_mcmc_iter_not$blockName;
arg_blocks_few_iterations <- arg_blocks_few_iterations[!arg_blocks_few_iterations %in% arg_block_names]



# subset by ARG block length ( >= min_block_len)
lines_mcmc_iter_long   <- lines_mcmc_iter[(lines_mcmc_iter$endInd - lines_mcmc_iter$startInd) >= min_block_len,]
lines_mcmc_iter_short  <- lines_mcmc_iter[(lines_mcmc_iter$endInd - lines_mcmc_iter$startInd) <  min_block_len,]
num_long_arg_blocks    <- nrow(lines_mcmc_iter_long)                  
num_short_arg_blocks   <- nrow(lines_mcmc_iter_short)                 
#                          total_length_contigs(lines_mcmc_iter)       
total_len_long_blocks  <- total_length_contigs(lines_mcmc_iter_long)  
total_len_short_blocks <- total_length_contigs(lines_mcmc_iter_short) 
arg_block_names        <- lines_mcmc_iter_long$blocknames;

# Write info into info file
write(paste("number of analyzed scaffolds is",num_scaffolds, sep= " " ), infoFile, append = FALSE )
write(paste("number of ARG blocks with at least",mcmc_iter,"MCMC iterations, of size at least",min_block_len," bp is",num_long_arg_blocks, sep= " " ), infoFile, append = TRUE )
write(paste("total length of these ARG blocks is:",total_len_long_blocks, sep= " " ), infoFile, append = TRUE )
write(paste("number of ARG blocks with at least",mcmc_iter,"MCMC iterations, but shorter than",min_block_len," bp is",num_short_arg_blocks, sep= " " ), infoFile, append = TRUE )
write(paste("total length of these ARG blocks is:",total_len_short_blocks, sep= " " ), infoFile, append = TRUE )

output=c()
for (scaffold in unique(lines_mcmc_iter_long$scaffold)){
  arg_blocks_per_scaffold <- lines_mcmc_iter_long[lines_mcmc_iter_long$scaffold == scaffold,]
  output = rbind(output , divideBlocksToRegionsBySegments(arg_blocks_per_scaffold, trim_len, segment_len));
}
write.table(output[,c("file","scaffold","startInd","endInd")], file = outFile, row.names= FALSE, quote = FALSE , col.names = FALSE, append = FALSE)






  ### R script to run on Linux  ###
############################################################################
##### compute average stats per window across entire saffold            ####
#####  - uses stats computes by tre_to_stat.R                           ####
#####  - averages stats across windows of given length (20 kb)          ####
#####  - produces BED file with average stats per window                ####
#####  - checks that all windows are covered equally by sampled trees   ####
#####  - writes information on scaffold coverage in end of given info file #
############################################################################

############################################################################
# script arguments and i/o
# example:   args         <- c("./argStats","Contig467", "20000", "argStats-windows", "infoTables/scaffold-window-coverage.txt")
statsDir     <- args[1]
scaffold     <- args[2]
window_len  <- as.numeric(args[3])
outDir       <- args[4]
infoFile     <- args[5]
# rounding parameter
ageDigits    <- 0   # keep integer value
enrichDigits <- 2   # keep 2 decimal digits in enrichment score
############################################################################

############################################################################
# create table with names of stats files for arg blocks in scaffold and their coordinates
statFiles     <- list.files(statsDir,pattern=paste(scaffold,"\\.",sep=""))
startPos      <- sub(paste(".*",scaffold,"\\.",sep=""), "", statFiles)
startPos      <- sub("-[0-9].*" ,"", startPos)
endPos        <- sub(".*-"     ,"", statFiles)
endPos        <- sub("\\..*"   ,"", endPos)

arg_block_table           <- data.frame(matrix(nrow=length(statFiles), ncol=3))
colnames(arg_block_table) <- c("file", "startPos", "endPos")
arg_block_table$file      <- paste(statsDir,statFiles,sep="/")
arg_block_table$startPos  <- as.numeric(startPos)
arg_block_table$endPos    <- as.numeric(endPos)
# order based on start coordinate of block
arg_block_table           <- arg_block_table[order(arg_block_table$startPos),]


############################################################################
# prepare output file according to columns of stat files
# (replacing first two columns with (chrom startPos endPos)
outFile        <- paste(outDir,"/",scaffold,".stat.bed",sep="")
dir.create(outDir,showWarnings=FALSE)
res <- file.create(outFile,force=TRUE)
treeStats      <- read.table(arg_block_table$file[1],quote="",header=TRUE)
stat_list      <- colnames(treeStats)[-(1:2)]
age_stats      <- c("TMRCA_all", "RT12", stat_list[grep("RTH$",stat_list)])
enrich_stats   <- stat_list[grep("_enrich$",stat_list)]
window_stats <- c("chrom", "startPos", "endPos", stat_list)
write.table(t(window_stats), file=outFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
# expected number of trees per window based on gap between first two stats
num_trees_per_seg <- window_len / (treeStats$pos[2] - treeStats$pos[1])
############################################################################

############################################################################
# main loop - iterate over arg block stat files
# pos holds current position and current_stretch_start holds start position
#   of current stretch of windows
pos <- arg_block_table$startPos[1]-1
current_stretch_start <- pos
for(i in 1:nrow(arg_block_table)) {
  if(i>1) {
    treeStats <- read.table(arg_block_table$file[i],quote="",header=TRUE)
  }
  chroms <- unique(treeStats$chrom)
  if(length(chroms) > 1 || chroms[1] != scaffold) {
     # unexpected name of chromosome in arg block
     message <- paste(scaffold,"found statistics associated with scaffolds",chroms, "in ARG block",i)
     write(message, file=infoFile, append=TRUE)
     break
  } else if(pos > arg_block_table$startPos[i]-1) {
     # unexpected overlap between arg blocks after trimming
     message <- paste(scaffold,"found overlap between ARG blocks",i-1,"( pos",pos,") and",i,"( pos",arg_block_table$startPos[i]-1,")")
     write(message, file=infoFile, append=TRUE)
     break
  } else if(pos < arg_block_table$startPos[i]-1) {
     # missing windows
     message <- paste(scaffold,current_stretch_start,pos,sep="\t")
     write(message,file=infoFile, append=TRUE)
     # cat("Found gap before arg block",i,"positions",pos,"-",arg_block_table$startPos[i]-1,"\n")
     pos <- arg_block_table$startPos[i]-1
     current_stretch_start <- pos
  }

  # loop over windows in arg block
  unexpected_halt <- FALSE
  while(pos < arg_block_table$endPos[i]) {
    window_lines <- which(treeStats$pos > pos & treeStats$pos <= pos+window_len)
    if(length(window_lines) != num_trees_per_seg) {
      # unexpected number of stat lines in window
      message <- paste(scaffold,"unexpected number of trees (",length(window_lines)," instead of",num_trees_per_seg,") in window strating at",pos," in ARG block",i)
      write(message, file=infoFile, append=TRUE)
      unexpected_halt <- TRUE
      break
    }
    # fill window stats - initialize chrom and positions
    window_stats_df <- data.frame(matrix(nrow=1,ncol=length(window_stats)))
    colnames(window_stats_df)    <- window_stats
    window_stats_df$chrom        <- scaffold
    window_stats_df$startPos     <- pos
    window_stats_df$endPos       <- pos+window_len
    # compute mean stats per trees in window and round off
    window_stats_df[1,stat_list] <- colMeans(treeStats[window_lines,stat_list])
    window_stats_df[1,enrich_stats] <- round(window_stats_df[1,enrich_stats] , digits=enrichDigits)
    window_stats_df[1,age_stats]    <- round(window_stats_df[1,age_stats]    , digits=ageDigits)
    
    # write stats line and advance pos
    write.table(window_stats_df[1,window_stats], file=outFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
    pos = pos + window_len
  } # end of while(pos)

  if(unexpected_halt) {
    break
  }

  if(pos > arg_block_table$endPos[i]) {
    # unexpected gap in the end
     message <- paste(scaffold," gap in end of ARG block",i,"( pos",arg_block_table$endPos[i],"-",pop,")")
     write(message, file=infoFile, append=TRUE)
     break
  }
} # end of for for(i)
############################################################################

############################################################################
# finalize - write final stretch in info file (should be one stretch per scaffold if all is okay)
message <- paste(scaffold,current_stretch_start,pos,sep="\t")
write(message,file=infoFile, append=TRUE)
# gzip file - TURNED OFF
# system(paste("gzip",outFile,"-f"))
############################################################################

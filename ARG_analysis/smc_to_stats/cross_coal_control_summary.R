# R script for summarizing results from cross coalescence statistics on control blocks


############################################################################
# set these
inFile   <- "infoTables/cross_coals_control_blocks.tsv"   # file with quantile differences summarized in all control blocks
outFile  <- "infoTables/cross_coals_control_blocks_summary.txt"  # file where text summary is written
pvalue   <- 0.01                                          # significance threthols
upperOrLowerTail <- TRUE                                  # compute upper tail (for CC elevation; TRUE) or lower tail (for CC reduction; FALSE)
############################################################################

ccTable     <- read.table(inFile,header=TRUE)
pairs       <- colnames(ccTable)
numRegions  <- nrow(ccTable)

tailLength  <- floor(numRegions*pvalue)
cat("thresholds for p-value = ",pvalue," using ",numRegions," control regions\n",file=control_summary_file, append=TRUE)
for(pair in pairs) {
  tail      <- sort(ccTable[,pair],decreasing=upperOrLowerTail)[1:(tailLength+1)]
  thres     <- mean(rev(tail)[1:2])
  cat(pair, ": (", thres, ") ", paste(sep=' ', tail),"\n")
}
sortedVals <- sort(as.matrix(ccTable), decreasing=upperOrLowerTail)
tailLength  <- floor(length(sortedVals)*pvalue)

cat("Threshold across pairs in p-value ",pvalue,": ",mean(sortedVals[tailLength+0:1])," extreme val: ", sortedVals[1],"\n")

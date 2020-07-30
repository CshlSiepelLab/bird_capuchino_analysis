# R script for summarizing results from cross coalescence statistics on control blocks


############################################################################
# set these
inFile   <- "crossCoals/cross_coals_control_v4.tsv"   # file with quantile differences summarized in all control blocks
pvalue   <- 0.01                                    # significance threthols
############################################################################

ccTable     <- read.table(inFile,header=TRUE)
pairs       <- colnames(ccTable)
numRegions  <- nrow(ccTable)

tailLength  <- floor(numRegions*pvalue)
cat("thresholds for p-value = ",pvalue," using ",numRegions," control regions\n")
for(pair in pairs) {
  tail      <- sort(ccTable[,pair],decreasing=TRUE)[1:(tailLength+1)]
  # tail      <- sort(ccTable[,pair],decreasing=FALSE)[1:(tailLength+1)]  # UNCOMMENT TO GET LEFT TAIL
  thres     <- mean(rev(tail)[1:2])
  cat(pair, ": (", thres, ") ", paste(sep=' ', tail),"\n")
}
sortedVals <- sort(as.matrix(ccTable), decreasing=TRUE)
# sortedVals <- sort(as.matrix(ccTable), decreasing=FALSE)  # UNCOMMENT TO GET LEFT TAIL
tailLength  <- floor(length(sortedVals)*pvalue)

cat("Threshold across pairs in p-value ",pvalue,": ",mean(sortedVals[tailLength+0:1])," extreme val: ", sortedVals[1],"\n")

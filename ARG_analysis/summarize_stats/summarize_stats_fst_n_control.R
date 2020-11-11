##########################################################################################################
###################    summarize max population enrichment and min RTH' in genomic block  ################
##########################################################################################################

############################################################################
# set these based on required empirical p-value
pvalue_stats                   <-  0.0001     # empirical p-value used for enrichment scores and RTH' (tests 1 & 2)
pvalue_RT12                    <-  0.001      # relaxed empirical p-value used for RT12
printingTail                   <-  0.0005     # fraction of most extreme values to print in report (should be > pvalue_stats)
# toggle TRUE/FALSE the three types of computations this script does
compute_FST_stats              <- TRUE        # anayze Fst peaks
compute_control_block_stats    <- TRUE        # analyze control blocks
compute_control_scaffold_stats <- TRUE        # analyze control windows
############################################################################


# input files
windowStatDir           <- "./argStats-windowed"
fst_peak_file           <- "./infoTables/FST-table-base.txt"
control_block_file      <- "./infoTables/control-blocks.txt"
control_scaffold_file   <- "./infoTables/control-scaffolds.txt"


# output files
fst_stats_file       <- "./infoTables/fst_peak_arg_stats.tsv"
control_stats_file   <- "./infoTables/control_blocks_arg_stats.tsv"
control_summary_file <- "./infoTables/control_regions_summary.txt"

# new source for F_ST peak table in shared directory
fst_peak_table <- read.table(fst_peak_file, header=TRUE, stringsAsFactors=FALSE)
# table with coordinates of all control blocks (scaffolds not containing F_ST peaks, split into 500 kb blocks)
if(compute_control_block_stats) {
   control_blocks <- read.table(control_block_file, stringsAsFactors=FALSE)
   colnames(control_blocks) <- c("scaffold", "startPos", "endPos")
}

# table with coordinates of all control scaffolds (scaffolds not containing F_ST peaks)
if(compute_control_scaffold_stats) {
   control_scaffolds <- read.table(control_scaffold_file, stringsAsFactors=FALSE)
   colnames(control_scaffolds) <- c("scaffold", "startPos", "endPos")
}

################################################
sumStats <- function(windowStatDir, scaffold, startPos, endPos){
  # outer boundaries are plotting boundaries and inner boundaries are for highlighted region (e.g. F_ST peak)
  
  # windowed stats file
  filename <- sprintf("%s/%s.stat.bed", windowStatDir, scaffold)
  stat_df <- as.data.frame(read.table(filename, header=TRUE, stringsAsFactors=FALSE))
  stats_all <- colnames(stat_df)
  # subset rows to take that intersect with the region
  df_peak_region <- stat_df[ stat_df$endPos >= startPos & stat_df$startPos <= endPos,] 

  # take a subset with the desired columns
  stats_enrich     <- grep("_enrich$",stats_all,value=TRUE)
  max_enrich       <- apply(df_peak_region[,stats_enrich] , 2, max)


  stats_RTH    <- grep("_RTH$",stats_all, value=TRUE)
  stats_RT     <- c(stats_RTH,"RT12")
  min_RT       <- round(digits=3, apply(df_peak_region[,stats_RT], 2, min))

  stats_sum                 <- c("scaffold", "startPos", "endPos", stats_enrich, stats_RT)
  block_stats_df            <- data.frame(matrix(nrow=1,ncol=length(stats_sum)))
  colnames(block_stats_df)  <- stats_sum
  block_stats_df[1,"scaffold"]  <- scaffold
  block_stats_df[1,"startPos"] <- startPos
  block_stats_df[1,"endPos"]   <- endPos
  block_stats_df[1,stats_enrich]       <- max_enrich
  block_stats_df[1,stats_RT]   <- min_RT

  return(block_stats_df)
}
################################################


################################################
# FST peaks
################################################
if(compute_FST_stats) {
   fst_stat_table <- c()
   for(i in 1:nrow(fst_peak_table)) {
      cat("summarizing stats for FST peak Scaffold",fst_peak_table[i,"peakID"],"\n")

      scaffold   <- paste("Contig",fst_peak_table[i,"Scaffold"], sep="")
      startPos   <- as.numeric(fst_peak_table[i,c("startPos")])
      endPos     <- as.numeric(fst_peak_table[i,c("endPos")])
      peak_stats <- sumStats(windowStatDir, scaffold, startPos, endPos)
      fst_stat_table <- rbind(fst_stat_table , peak_stats)
   } # end of for(i)
   rownames(fst_stat_table) <- fst_peak_table[,"peakID"]
   write.table(fst_stat_table, fst_stats_file, quote=FALSE, sep="\t")
}


################################################
# control blocks
################################################
if(compute_control_block_stats) {
   cat("summarizing stats for control regions (. per 50 regions)")

   control_stat_table <- c()
   for(i in 1:nrow(control_blocks)) {
      scaffold   <- control_blocks[i,"scaffold"]
      startPos   <- as.numeric(control_blocks[i,c("startPos")])
      endPos     <- as.numeric(control_blocks[i,c("endPos")])
      peak_stats <- sumStats(windowStatDir, scaffold, startPos, endPos)
      control_stat_table <- rbind(control_stat_table , peak_stats)
      if(i %% 50 == 0) cat(".")
   } # end of for(i)
   cat("\n")
   write.table(control_stat_table, control_stats_file, quote=FALSE, sep="\t", row.names=FALSE)
}


################################################
# control windows
################################################
if(compute_control_scaffold_stats) {
   cat("collecting stats for control scaffolds (. per 50 scaffold)")

   control_scaffold_table <- c()
   for(i in 1:nrow(control_scaffolds)) {
      scaffold   <- control_scaffolds[i,"scaffold"]
      startPos   <- as.numeric(control_scaffolds[i,c("startPos")])
      endPos     <- as.numeric(control_scaffolds[i,c("endPos")])

      filename <- sprintf("%s/%s.stat.bed", windowStatDir, scaffold)
      stat_df <- as.data.frame(read.table(filename, header=TRUE, stringsAsFactors=FALSE))
      # subset rows to take that intersect with the region
      df_peak_region <- stat_df[ stat_df$endPos >= startPos & stat_df$startPos <= endPos,]

      control_scaffold_table <- rbind(control_scaffold_table , df_peak_region)
      if(i %% 50 == 0) cat(".")
   } # end of for(i)
   cat("\n")

   ################################################
   # Generate report:
   ################################################

   # show most extreme values for every stat
   numVals   <- round(nrow(control_scaffold_table) * printingTail, digits=0)
   ind_pval  <- round(nrow(control_scaffold_table) * pvalue_stats - 0.5, digits=0)
   stats_enrich     <- grep("_enrich$",colnames(stat_df),value=TRUE)
   stats_RTH <- grep("_RTH$",colnames(stat_df),value=TRUE)
   stats_RT  <- c(stats_RTH,"RT12")
   cat("Summary for ",nrow(control_scaffold_table)," control windows in non-FST scaffolds:\n",file=control_summary_file)
   thresholds <- c()
   for(stat in stats_enrich) {
      vals <- sort(control_scaffold_table[,stat],decreasing=TRUE)[1:numVals]
      cat("Top ",numVals," values for ",stat," (out of ",nrow(control_scaffold_table),") are ", vals,"\n",file=control_summary_file, append=TRUE)
      # setting individual threhold per pop for enrichment
      threshold <- round(mean(vals[0:1+ind_pval]),digits=1)
      thresholds[stat] <- threshold
      # cat("  --> Setting threshold for ",stat," to ",threshold,"\n",file=control_summary_file, append=TRUE)
   }
   for(stat in stats_RTH) {
      vals <- sort(control_scaffold_table[,stat],decreasing=FALSE)[1:numVals]
      cat("Bottom ",numVals," values for ",stat," (out of ",nrow(control_scaffold_table),") are ", vals,"\n",file=control_summary_file, append=TRUE)
      # determine species-specific thresholds
      threshold <- round(mean(vals[0:1+ind_pval]),digits=1)
      thresholds[stat] <- threshold
      # cat("  --> Setting threshold for ",stat," to ",threshold,"\n",file=control_summary_file, append=TRUE)
   }

   # separate p-value threshold for RT12
   stat      <- "RT12"
   ind_RT12 <- round(nrow(control_scaffold_table) * pvalue_RT12 - 0.5, digits=0)
   numVals   <- ind_RT12+10
   vals <- sort(control_scaffold_table[,stat],decreasing=FALSE)[1:numVals]
   cat("Bottom ",numVals," values for ",stat," (out of ",nrow(control_scaffold_table),") are ", vals,"\n",file=control_summary_file, append=TRUE)
   # setting individual threhold per pop for enrichment
   threshold <- round(mean(vals[0:1+ind_RT12]),digits=3)
   thresholds[stat] <- threshold
   # cat("  --> Setting threshold for ",stat," to ",threshold,"\n",file=control_summary_file, append=TRUE)

   cat("\nThresholds:\n",file=control_summary_file, append=TRUE)
   write.table(thresholds,quote=FALSE,col.names=FALSE,sep="\t",file=control_summary_file, append=TRUE) 

   if(compute_control_block_stats) {   
      cat("\nChecking threhsolds on control blocks:\n",file=control_summary_file, append=TRUE)
      for(stat in stats_enrich) {
         count <- length(which(control_stat_table[,stat] >= thresholds[stat]))
         cat("There are ", count, " control blocks with ",stat," > ",thresholds[stat]," (out of ", nrow(control_stat_table),")\n" ,file=control_summary_file, append=TRUE)
      }
      for(stat in stats_RT) {
         count <- length(which(control_stat_table[,stat] <= thresholds[stat]))
         cat("There are ", count, " control blocks with ",stat," < ",thresholds[stat]," (out of ", nrow(control_stat_table),")\n" ,file=control_summary_file, append=TRUE)
      }
   }
}
warnings()


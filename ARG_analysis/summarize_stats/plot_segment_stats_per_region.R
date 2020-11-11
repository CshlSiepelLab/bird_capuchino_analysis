#########################################################################################################
############  plot population enrichment, RTH, and cross coalescence along a scaffold region  ###########
############  uses prepared windowed-stat files for this                                      ###########
#########################################################################################################

# This script plots the enrichment scores RTH values and normalized cross coalescence times onto a pdf file
# inputs: a path to a dir with windowed statistics and a coordinates of region to plot.
#         flags can be used to enable / disable some plots
# output: A pdf file with all plots. 
# The order of plots: enrichment scores, RTH, cross coalescence.


library(reshape2) # for melt function
library(ggplot2)
library(cowplot)


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
# print(script_dir)

############################################################################
# args <- c("./argStats-windowed" , "./plots/testPlot.pdf" , "Contig252" , "220000" , "710000", "420000" , "510000") #EXAMPLE ARGUMENT LINE FOR PEAK ON Contig252
args     <- commandArgs(trailingOnly = TRUE) # TRUE
segmentStatDir <- args[1]
plotFile       <- args[2]
scaffold       <- args[3]
startPos       <- as.numeric(args[4])
endPos         <- as.numeric(args[5])
startRegion    <- as.numeric(args[6])      # region highlighted with vertical dashed bars (Fst peak) (set to -1 for no region highlight)
endRegion      <- as.numeric(args[7])      # region highlighted with vertical dashed bars (Fst peak) (set to -1 for no region highlight)
hasRegion      <- FALSE
if(startRegion >= startPos && endRegion <= endPos) {
   hasRegion <- TRUE
}
############################################################################


############################################################################
# script setup
# flags for plotting
plot_enrich     <- TRUE
plot_rth        <- TRUE
plot_crosscoal  <- TRUE
# threshold for species enrichment based on value that best represents the significance threshold across species.
enrichThres <- 4.0
rthThres <- 0.4
crosscoalThres <- -1

# plot height per panel (in inches) and width proportional to range length
maxRangeLength <- 1000000
panelHeight     <- 3
plotWidth    <- panelHeight*5*(endPos-startPos)/maxRangeLength

# pop_ind_key used to assign stats with colors
pop_ind_key   <- read.csv("infoTables/individual-species-key-modified.txt",sep="\t", stringsAsFactors = FALSE)
############################################################################

pop_col_table <- unique(pop_ind_key[,c("Species","Color")])
populations   <- pop_col_table[,"Species"]
# use dark gray for the RT12 stat
colors        <- c("darkgray",pop_col_table[,"Color"])
names(colors) <- c("RT12",populations)

# get the bed stat file to be read with RTH and species enrichment stats
# and filter in only rows of stat file that intersect the region
filename <- sprintf("%s/%s.stat.bed", segmentStatDir, scaffold)
stat_df <- as.data.frame(read.table(filename, header=TRUE, stringsAsFactors=FALSE))
stat_df <- stat_df[ stat_df$endPos > startPos & stat_df$startPos < endPos,] 

# do the same with all stat files for cross coalescence for the relevant scaffold
filenames <- paste(sep="/",segmentStatDir,list.files(path=segmentStatDir, pattern=paste(scaffold,"[.].*[.]crosscoal.bed",sep="")))
crosscoal_df <- c()
for(file in filenames) {
   crosscoal_df <- rbind(crosscoal_df ,  as.data.frame(read.table(file, header=TRUE, stringsAsFactors=FALSE)))
}
crosscoal_df <- crosscoal_df[ crosscoal_df$endPos > startPos & crosscoal_df$startPos < endPos,] 

plotHeight <- 0
plot_list <- list()
plot_heights <- c()
if(plot_enrich) {
  cat("Generating species enrichment plot.\n")
  
  plotHeight  <- plotHeight + panelHeight
  ##########################################################
  # plopulation enrichment plot
  ##########################################################
  # take a subset with the desired columns
  enrich_df <- stat_df[,c("startPos",paste(populations,"enrich",sep="_"))]
  colnames(enrich_df) <-   c("startPos", populations)
  # melt the df by startPos  
  melted_stat_df <- melt(enrich_df, id = "startPos") #  

  options(repr.plot.width = plotWidth,repr.plot.height = panelHeight)
  enrichment_plot <- ggplot(melted_stat_df, aes(x = as.numeric(startPos), y =as.numeric(value), group =variable , colour = variable))
  enrichment_plot <- enrichment_plot + geom_point( alpha = 0.75)
    # colors and legend
  enrichment_plot <- enrichment_plot + scale_colour_manual(values = colors[populations])
  enrichment_plot <- enrichment_plot + labs(color='')
    # remove grid and background
  enrichment_plot <- enrichment_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
    # axes
  enrichment_plot <- enrichment_plot + scale_y_continuous( expand = c(0, 0),name="species\nenrichment", limits = c(1, 11) )
  enrichment_plot <- enrichment_plot + scale_x_continuous(expand = c(0, 0),name="position on scaffold" , limits = c(startPos,endPos), labels = scales::comma)
  if(plot_rth || plot_crosscoal) {
     enrichment_plot <- enrichment_plot + theme(axis.title.x=element_blank()) #,  axis.text.x=element_blank(), axis.ticks.x=element_blank())
  } else {
     enrichment_plot <- enrichment_plot + theme(axis.text.x = element_text(size=6, angle=0))
  }
  # legend
  if(plot_rth) {
     enrichment_plot <- enrichment_plot + theme(legend.position="none")  
     plot_heights              <- c(plot_heights , 1)
  } else {
     enrichment_plot <- enrichment_plot + theme(legend.position="bottom")
     plot_heights              <- c(plot_heights , 1.3)
  }
    ## add horizontal lines at thresholds
  if(enrichThres != -1) { 
     enrichment_plot <- enrichment_plot + geom_hline(yintercept = enrichThres, linetype = "dashed", color = "darkgray", size = 0.3)
  }
  
  # add  lines for peak range
  if(hasRegion) {
     enrichment_plot <- enrichment_plot + geom_vline(xintercept = c(startRegion,endRegion), linetype="dotted", color = "blue")
  }
  plot_list$enrichment_plot <- enrichment_plot
}# end of if(plot_enrich)

if(plot_rth) {
  cat("Generating RTH plot.\n")
  plotHeight  <- plotHeight + panelHeight

  ##########################################################
  # RTH plot
  ##########################################################
  # take a subset with the desired columns and normalize clade12 ages
  rth_df <- stat_df[,c("startPos",paste(populations,"RTH",sep="_"), "RT12")]
  # sub_df_peak_region[,2:ncol(sub_df_peak_region)] <- sub_df_peak_region[,2:ncol(sub_df_peak_region)] 
  colnames(rth_df) <-   c("startPos", populations, "RT12")
  # melt the df by the startind coordinate column  
  melted_stat_df <- melt(rth_df, id = "startPos") #  

  options(repr.plot.height = plotWidth,repr.plot.height = panelHeight)
  rth_plot <- ggplot(melted_stat_df, aes(x = as.numeric(startPos), y =as.numeric(value), group =variable , colour = variable)) 
  rth_plot <- rth_plot + geom_point( alpha = 0.75)
    # colors and legend (none)
  rth_plot <- rth_plot + scale_colour_manual(values = c(colors[populations],"darkgray"))
  rth_plot <- rth_plot + labs(color='')
    # remove grid and background
  rth_plot <- rth_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
    # axes
  rth_plot <- rth_plot + scale_y_continuous( expand = c(0, 0),name="RTH", limits = c(0.0,1.1) )
  rth_plot <- rth_plot + scale_x_continuous(expand = c(0, 0),name="Position on scaffold" , limits = c(startPos,endPos), labels = scales::comma)
  if(plot_crosscoal) {
     rth_plot <- rth_plot + theme(axis.title.x=element_blank()) #, axis.text.x=element_blank(), axis.ticks.x=element_blank())
  } else {
     rth_plot <- rth_plot + theme(axis.text.x = element_text(size=6, angle=0))
  }
  # legend
  rth_plot <- rth_plot + theme(legend.position="bottom")
  
    ## add horizontal lines at thresholds
  if(rthThres != -1) {
     rth_plot <- rth_plot + geom_hline(yintercept = rthThres,linetype = "dashed", color = "darkgray", size = 0.3)
  }

  # add  lines for peak range
  if(hasRegion) {
     rth_plot <- rth_plot + geom_vline(xintercept = c(startRegion,endRegion), linetype="dotted", color = "blue")
  }
  plot_list$rth_plot <- rth_plot
  plot_heights              <- c(plot_heights , 1.3)
}# end of if(plot_rth)


if(plot_crosscoal) {  
  cat("Generating cross coalescence plot.\n")
  plotHeight  <- plotHeight + panelHeight

  ##########################################################
  # Cross coalescence plot
  ##########################################################
  # population pairs
  popPairs     <- c()
  for(i in 1:(length(populations)-1)) {
    for(j in (i+1):length(populations)) {
      popPairs <- c(popPairs , paste(populations[i],populations[j],sep="_"))
    }
  }
  # discard the contig name column because it is a different type of attribute
  melted_stat_df <- melt(crosscoal_df[,c("startPos",popPairs)] , id ="startPos") 
  
  options(repr.plot.height = plotWidth,repr.plot.height = panelHeight)
  crosscoal_plot <- ggplot(melted_stat_df , aes(x = as.numeric(startPos), y =as.numeric(value), group =variable , colour = variable))
  crosscoal_plot <- crosscoal_plot + geom_point( alpha = 1.0) + scale_colour_manual(values = rainbow(length(popPairs)))
  crosscoal_plot <- crosscoal_plot + ylab("Normalized cross\ncoalescence times")
  # remove grid and background
  crosscoal_plot <- crosscoal_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
  crosscoal_plot <- crosscoal_plot  + theme(axis.text.x = element_text(size=6, angle=0)) + theme(legend.position="bottom", legend.title = element_blank())
  crosscoal_plot <- crosscoal_plot  + guides(colour=guide_legend(ncol=4))
  crosscoal_plot <- crosscoal_plot + scale_x_continuous(expand = c(0, 0),name="Position on scaffold" , limits = c(startPos,endPos), labels = scales::comma)
  # add  lines for peak range
  if(hasRegion) {
     crosscoal_plot <- crosscoal_plot + geom_vline(xintercept = c(startRegion,endRegion), linetype="dotted", color = "blue")
  }
  plot_list$crosscoal_plot <- crosscoal_plot
  plot_heights              <- c(plot_heights , 1.5)
}# end of if(plot_crosscoal)


if(plotHeight > 0) {
   cat("Creating plot file ",plotFile,"\n")
   pdf(file = plotFile ,width = plotWidth, height = plotHeight)
   print(cowplot::plot_grid(plotlist=plot_list, ncol=1, align="v" , axis = "l" ,rel_heights = plot_heights) , scale = c(1, .9, .9, .7))
   dev.off()
} else {
   cat("No plots plotted. Set plotting flags to enable plotting.\n")
}


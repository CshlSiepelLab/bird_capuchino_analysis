library(ape)
library(phytools)
library(ggtree)

args <- commandArgs(TRUE)

id_ <- as.integer(args[1]) #scaffold id
path_to_species_ids <- args[2] #file containing species ids
segment_num <- args[3] #segment run
path_input = args[4] #path to input file
path_output = args[5] #path to output file
path_positions = args[6] #path to file containing positions to plot: start, end

data <- read.table(path_to_species_ids)
id <- c()
for(i in 1:dim(data)[1]){
  id <- rbind(id, cbind(toString(data[i,1]), paste(toString(data[i,2]), "_", "1", sep="")))
  id <- rbind(id, cbind(toString(data[i,1]), paste(toString(data[i,2]), "_", "2", sep="")))
}

table <- read.table(paste(path_input, "/Contig", id_,"-trees/out.1000_", segment_num,".bed.gz", sep=""))
table <- table[,2:3]
tree1 <- read.tree(paste(path_input, "/Contig", id_,"-trees/out.1000_", segment_num,".bed.gz", sep=""))

table_pos <- read.table(path_positions)
for(kk in 1:dim(table_pos)[1]){
  min_ <- table_pos[kk,1]*1e6                                                                    
  max_ <- table_pos[kk,2]*1e6                                                                    
  table[which(table[,1]>=min_ & table[,2]<=max_),]
  index_ <- which(table[,1]>=min_ & table[,2]<=max_)
  POS_ <- table[which(table[,1]>=min_ & table[,2]<=max_),]
  for(iter in 1:length(index_)){
    j <- index_[iter]
    tree <- tree1[[j]]
    val <- c()
    for(i in 1:length(tree$tip.label)){
      index <- which(id[,2]==tree$tip.label[i])
      if(length(index)==0){
        val <- c(val, tree$tip.label[i])
        #tree$tip.label[i] <- ""
      } else {
        tree$tip.label[i] <- id[index, 1]
      }
    }
    
    for(i in val){
      print(i)
      tree <- drop.tip(tree, i)
    }
    
    v <- c()
    code <- c()
    for(k in tree$tip.label){
      if(length(grep('mel', k))>0){
        code <- c(code, "mel")
      } else if(length(grep('nig', k))>0){
        code <- c(code, "nig")
      } else if(length(grep('hyp', k))>0){
        code <- c(code, "hyp")
      } else if(length(grep('pal', k))>0){
        code <- c(code, "pal")
      } else if(length(grep('pil', k))>0){
        code <- c(code, "pil")
      }
    }
    
    dd <- data.frame(tree$tip.label, code)
    row.names(dd) <- NULL
    colnames(dd) <- c("taxa", "Species")
    print(dd)
    plot <- ggtree(tree) %<+% dd + geom_tippoint(aes(color=Species)) + theme(legend.position="right", text = element_text(size=26) ) + scale_color_manual(values=c("red", "blue", "green", "magenta", "orange")) #+ geom_hilight(node=228, fill="green") + geom_hilight(node=161, fill="blue") + guides(colour = guide_legend(override.aes = list(size=5))) + geom_cladelabel(node=161, label="", color="blue", offset=.2, align=TRUE) + geom_cladelabel(node=228, label="", color="green", offset=.2, align=TRUE) 

    min_ <- POS_[iter,1]
    max_ <- POS_[iter,2]
    filename=paste(path_output, "/Contig", id_,"-trees/tree_", id_, "_", min_, "_", max_, "_", iter, ".pdf", sep="")
    pdf(file = filename)
    par(mar=c(4.5,5,1,1))
    print(plot)
    dev.off()
  }
}

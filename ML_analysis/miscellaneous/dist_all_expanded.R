library(lattice)
library(qqman)
library("CMplot")
path <- "expanded"

setwd(paste('~/Desktop/', path, sep=""))

seg <- read.table("seg.txt")
seg <- seg[[1]]
for(i in 2:length(seg)){
  seg[i] <- seg[i] + seg[i-1]
}

maxY <- 1
ids = c(118,1635,1717,1954,252,257,263,3622,404,412,430,567,579,59,637,762,766,791)
type <- c("hypox", "mel", "nig", "pal", "pil")

contigs <- c()
for(folder in type){
  start <- 1
  counter <- 1
  dd <- c()
  
  for(id in ids){
    threshold <- 1000000
    contig <- ""
    file = paste("~/Desktop/", path, "/empirical_soft_test.txt",sep="")
    data <- read.table(file)
    
    end <- seg[counter]
    
    file = paste("~/Desktop/", path, "/pos.txt",sep="")
    pos <- read.table(file)
    pos <- pos[start:end,]
    c1 <- (pos[,2]-pos[,1])*0.5 + pos[,1]
    
    s1 <- data[start:end,2]
    s2 <- data[start:end,3]
    s3 <- data[start:end,4]
    s4 <- data[start:end,5]
    s5 <- data[start:end,6]
    
    if(folder==type[1])
      pval <- s1
    if(folder==type[2])
      pval <- s2
    if(folder==type[3])
      pval <- s3
    if(folder==type[4])
      pval <- s4
    if(folder==type[5])
      pval <- s5
    
    pval[which(pval==1)] <- 0.99
    dd <- rbind(dd, data.frame(c1, rep.int(id, length(pval)), c1, pval, rep.int(length(pval), length(pval))))
    
    start <- seg[counter]+1
    counter <- counter + 1
  }
  
  colnames(dd) <- c("SNP", "Chromosome", "Position", "P", "size")
  dd <- dd[order(dd$size),]
  dd$Chromosome <- dd$size
  dd <- dd[,1:4]
  dd$Position <- dd$Position/threshold
  CMplot(dd, plot.type="m",band=0, LOG10=FALSE, ylab="Prob", xlab = "Scaffolds", threshold.lty=1, threshold.lwd=1, threshold.col="red", amplify=FALSE, signal.col=NULL, chr.den.col=NULL, file="jpg", memo=folder, dpi=300)
}
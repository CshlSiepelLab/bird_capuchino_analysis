library(lattice)

type1 <- args[1] #species number 1 e.g. mel
type2 <- args[2] #species number 2 e.g. nig
path <- args[3] #species pair e.g. mel-nig
dir <- 'hardsoft'

setwd(paste('~/Desktop/', path, sep=""))

seg <- read.table("seg.txt")
seg <- seg[[1]]
for(i in 2:length(seg)){
  seg[i] <- seg[i] + seg[i-1]
}
sum_ <- 0
maxY <- 1

ids = c(118,1635,1717,1954,252,257,263,3622,404,412,430,567,579,59,637,762,766,791)
type <- c(paste("Complete soft (", type1,")", sep=""), paste("Complete soft (", type2,")", sep=""), paste("Complete hard (", type1,")", sep=""), paste("Complete hard (", type2,")", sep=""))
start <- 1
counter <- 1
for(id in ids){
  jpeg(file = paste("~/Desktop/", path, "/", dir, "/Contig", id, ".jpg", sep=""))
  par(mar=c(6,6.5,2,2))
  for(folder in type){
    threshold <- 1000000
    contig <- ""
    file = paste(dir, "/empirical_", dir, "_", path, ".txt", sep="")
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
    
    if(folder==type[1])
      pval <- s1
    if(folder==type[2])
      pval <- s2
    if(folder==type[3])
      pval <- s3
    if(folder==type[4])
      pval <- s4
    
    dd <- data.frame(rep.int(id, length(pval)), c1, pval)
    colnames(dd) <- c("chr", "pos", "pvalue")
    dd$pos <- dd$pos/threshold
    
    if(folder==type[1]){
      plot(dd$pos, dd$pvalue, pch=20, cex.axis=1.5, cex.main=2, cex.lab=2, xlim = c(2.5,6.5), ylim = c(0, maxY), col="brown", xlab = "Genomic position (Mb)", main = paste("Scaffold ", id, sep=""), ylab="Prob", type="n")
      Y <- c(0, 1)
      if(id==252){
        X <- c(0.42, 0.5)
        legend("topright", cex=1.25, legend = c(type[1], type[2], type[3], type[4]), pch=16, col = c("blue", "green", "brown", "magenta"))
        rect(X[1], Y[1], X[2], Y[2], border = "white", col = "grey")
      } else if(id==412){
        X <- c(3.4, 3.55)
        X <- c(3.2, 3.75)
        rect(X[1], Y[1], X[2], Y[2], border = "white", col = "grey")
      } else if(id==404){
        X <- c(5.05, 5.815)
        rect(X[1], Y[1], X[2], Y[2], border = "white", col = "grey")
      } else if(id==257){
        X <- c(21.62, 21.72)
        X1 <- c(24.4, 24.5)
        X <- c(21.42, 21.92)
        X1 <- c(24.2, 24.7)
        rect(X[1], Y[1], X[2], Y[2], border = "white", col = "grey")
        rect(X1[1], Y[1], X1[2], Y[2], border = "white", col = "grey")
      } else if(id==1717){
        X <- c(0.94, 1)
        rect(X[1], Y[1], X[2], Y[2], border = "white", col = "grey")
      } else if(id==263){
        X <- c(0.05, 0.35)
        rect(X[1], Y[1], X[2], Y[2], border = "white", col = "grey")
      } else if(id==791){
        X <- c(8.5, 10)
        rect(X[1], Y[1], X[2], Y[2], border = "white", col = "grey")
      }
      lines(dd$pos, dd$pvalue, pch=20, cex=1.2, col="blue", type = "p")
    } else if(folder==type[2]) {
      lines(dd$pos, dd$pvalue, pch=20, cex=1.2, col="green", type = "p")
    } else if(folder==type[3]) {
      lines(dd$pos, dd$pvalue, pch=20, cex=1.2, col="brown", type = "p")
    } else if(folder==type[4]) {
      lines(dd$pos, dd$pvalue, pch=20, cex=1.2, col="magenta", type = "p")
    }
  }
  
  start <- seg[counter]+1
  counter <- counter + 1
  dev.off()
}


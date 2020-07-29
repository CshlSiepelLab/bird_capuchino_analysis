library(pROC)
library(PRROC)

compute_ROC <- function(x,y,index){
  outcome <- c(rep.int(0,dim(x)[1]), rep.int(1,dim(y)[1]))
  data <- rbind(x, y)
  data <- data[,index]
  Test <- data.frame(outcome, data)
  colnames(Test) <- c("outcome", "prob")
  ROC1 <- roc(Test$outcome, Test$prob)
  print(ROC1$auc)
  return(ROC1)
}

type <- args[1] #species pair e.g. mel-nig
sp1 <- args[2] #species number 1 e.g. mel
sp2 <- args[3] #species number 2 e.g. nig
dir <- args[4] #type e.g. ancestral or soft

setwd(paste("~/Desktop/", type, sep=""))

if(dir=="parallel"){
  file = paste(dir, "/proba_parallel", type, ".txt",sep="")
  annot <- c(paste("neutral - soft (", sp1, ")", sep=""), paste("neutral - soft (", sp2, ")", sep=""), paste("neutral - parallel soft", sep=""), paste("soft (", sp1, ") - soft (", sp2, ")", sep=""), paste("soft (", sp1, ") - parallel soft", sep=""), paste("soft (", sp2, ") - parallel soft", sep=""))
}

if(dir=="ancestral"){
  sp1 <- "#1"
  sp2 <- "#2"
  file = paste(dir, "/proba_ancestral", type, ".txt",sep="")
  annot <- c(paste("neutral - soft (", sp1, ")", sep=""), paste("neutral - soft (", sp2, ")", sep=""), paste("neutral - ancestral soft", sep=""), paste("soft (", sp1, ") - soft (", sp2, ")", sep=""), paste("soft (", sp1, ") - ancestral soft", sep=""), paste("soft (", sp2, ") - ancestral soft", sep=""))
}

data <- read.table(file)

x1 <- 1000
x2 <- 1000
x3 <- 1000
x4 <- 1000 

neutral <- data[1:x1,]
soft1 <- data[(x1+1):(x1+x2),]
soft2 <- data[(x1+x2+1):(x1+x2+x3),]
ancestral <- data[(x1+x2+x3+1):(x1+x2+x3+x4),]

pdf(file = paste(dir, "/ROC_", type ,".pdf", sep=""))
par(mar=c(4, 4.5, 2.5, 0.5))
ROC1 <- compute_ROC(neutral, soft1, 2)
plot(1-ROC1$specificities, ROC1$sensitivities, cex.lab=2, cex.axis=2, cex.main=2, col="blue", xlab = "FPR", ylab = "TPR", type="n", pch=20, cex=2)
lines(1-ROC1$specificities, ROC1$sensitivities, col="antiquewhite", type="l", lwd=2)
ROC1 <- compute_ROC(neutral, soft2, 3)
lines(1-ROC1$specificities, ROC1$sensitivities, col="aquamarine", type="l", lwd=2)
ROC1 <- compute_ROC(neutral, ancestral, 4)
lines(1-ROC1$specificities, ROC1$sensitivities, col="aquamarine4", type="l", lwd=2)
ROC1 <- compute_ROC(soft1, soft2, 2)
lines(1-ROC1$specificities, ROC1$sensitivities, col="azure4", type="l", lwd=2)
ROC1 <- compute_ROC(soft1, ancestral, 2)
lines(1-ROC1$specificities, ROC1$sensitivities, col="black", type="l", lwd=2)
ROC1 <- compute_ROC(soft2, ancestral, 3)
lines(1-ROC1$specificities, ROC1$sensitivities, col="blue", type="l", lwd=2)
legend(0.59, 0.95, cex = 1, pch = 16, legend = c(annot[1], annot[2], annot[3], annot[4], annot[5], annot[6]), col = c("antiquewhite", "aquamarine", "aquamarine4", "azure4", "black", "blue"))
dev.off()

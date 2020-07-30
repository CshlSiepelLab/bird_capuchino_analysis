library(pROC)
library(PRROC)

args <- commandArgs(TRUE)

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
dir <- args[4] #type e.g. partial or hardsoft

setwd(paste("~/Desktop/", type, sep=""))

if(dir=="hardsoft"){
  file = paste(dir, "/proba_hardsoft", type, ".txt",sep="")
  sp1 <- "#1"
  sp2 <- "#2"
  annot <- c(paste("neutral - complete soft (", sp1, ")", sep=""), paste("neutral - complete soft (", sp2, ")", sep=""), paste("neutral - complete hard (", sp1, ")", sep=""), paste("neutral - complete hard (", sp2, ")", sep=""), paste("complete soft (", sp1, ") - complete soft (", sp2, ")", sep=""), paste("complete soft (", sp1, ") - complete hard (", sp1, ")", sep=""), paste("complete soft (", sp1, ") - complete hard (", sp2, ")", sep=""), paste("complete soft (", sp2, ") - complete hard (", sp1, ")", sep=""), paste("complete soft (", sp2, ") - complete hard (", sp2, ")", sep=""), paste("complete hard (", sp1, ") - complete hard (", sp2, ")", sep=""))
}

if(dir=="partial"){
  file = paste(dir, "/proba_partial", type, ".txt",sep="")
  sp1 <- "#1"
  sp2 <- "#2"
  annot <- c(paste("neutral - complete soft (", sp1, ")", sep=""), paste("neutral - complete soft (", sp2, ")", sep=""), paste("neutral - partial soft (", sp1, ")", sep=""), paste("neutral - partial soft (", sp2, ")", sep=""), paste("complete soft (", sp1, ") - complete soft (", sp2, ")", sep=""), paste("complete soft (", sp1, ") - partial soft (", sp1, ")", sep=""), paste("complete soft (", sp1, ") - partial soft (", sp2, ")", sep=""), paste("complete soft (", sp2, ") - partial soft (", sp1, ")", sep=""), paste("complete soft (", sp2, ") - partial soft (", sp2, ")", sep=""), paste("partial soft (", sp1, ") - partial soft (", sp2, ")", sep=""))
}  

data <- read.table(file)

x1 <- 1000
x2 <- 1000
x3 <- 1000
x4 <- 1000 
x5 <- 1000

neutral <- data[1:x1,]
soft1 <- data[(x1+1):(x1+x2),]
soft2 <- data[(x1+x2+1):(x1+x2+x3),]
hard1 <- data[(x1+x2+x3+1):(x1+x2+x3+x4),]
hard2 <- data[(x1+x2+x3+x4+1):(x1+x2+x3+x4+x5),]

pdf(file = paste(dir, "/ROC_", type ,".pdf", sep=""))
par(mar=c(4, 4.5, 2.5, 0.5))
ROC1 <- compute_ROC(neutral, soft1, 2)
plot(1-ROC1$specificities, ROC1$sensitivities, cex.lab=2, cex.axis=2, cex.main=2, col="blue", xlab = "FPR", ylab = "TPR", type="n", pch=20, cex=2)
lines(1-ROC1$specificities, ROC1$sensitivities, col="antiquewhite", type="l", lwd=2)
ROC1 <- compute_ROC(neutral, soft2, 3)
lines(1-ROC1$specificities, ROC1$sensitivities, col="aquamarine", type="l", lwd=2)
ROC1 <- compute_ROC(neutral, hard1, 4)
lines(1-ROC1$specificities, ROC1$sensitivities, col="aquamarine4", type="l", lwd=2)
ROC1 <- compute_ROC(neutral, hard2, 5)
lines(1-ROC1$specificities, ROC1$sensitivities, col="azure4", type="l", lwd=2)
ROC1 <- compute_ROC(soft1, soft2, 2)
lines(1-ROC1$specificities, ROC1$sensitivities, col="black", type="l", lwd=2)
ROC1 <- compute_ROC(soft1, hard1, 2)
lines(1-ROC1$specificities, ROC1$sensitivities, col="blue", type="l", lwd=2)
ROC1 <- compute_ROC(soft1, hard2, 2)
lines(1-ROC1$specificities, ROC1$sensitivities, col="brown", type="l", lwd=2)
ROC1 <- compute_ROC(soft2, hard1, 3)
lines(1-ROC1$specificities, ROC1$sensitivities, col="chartreuse", type="l", lwd=2)
ROC1 <- compute_ROC(soft2, hard2, 3)
lines(1-ROC1$specificities, ROC1$sensitivities, col="darkgoldenrod1", type="l", lwd=2)
ROC1 <- compute_ROC(hard1, hard2, 4)
lines(1-ROC1$specificities, ROC1$sensitivities, col="firebrick1", type="l", lwd=2)
legend(0.4, 0.75, cex = 1, pch = 16, legend = c(annot[1], annot[2], annot[3], annot[4], annot[5], annot[6], annot[7], annot[8], annot[9], annot[10]), col = c("antiquewhite", "aquamarine", "aquamarine4", "azure4", "black", "blue", "brown", "chartreuse", "darkgoldenrod1", "firebrick1"))
dev.off()

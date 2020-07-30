args <- commandArgs(TRUE)

"
Example run: Rscript find_rth_max.R pil 412 3501049 3511049 <path>
"

species <- args[1] #species name
id_1 <- as.integer(args[2]) #scaffold id
a <- as.integer(args[3]) #start position
b <- as.integer(args[4]) #end position		
path = args[5] #path to file
TEMP <- c()
arg <- c()
for(iter in 0:32){
  arg <- rbind(arg, read.table(paste(path, "/Contig", id_1,"." ,iter, ".stats.", species, ".txt", sep="")))
  TEMP = arg[which(arg[, 2] >= a & arg[, 3] <= b),]
  print(dim(TEMP))
  if(dim(TEMP)[1]>0){
    print(TEMP[min(TEMP[,7])==TEMP[,7],])
    break;
  }
}

####  R script to run in Linux  ###
#########################################################################################################
## read directory with smc files and log file from a log files directory  ##
## creates a new directory (if not exists) and returns bed files in the directory ##
## create .tbi files  in the same directory  ##  
## Command line:
### bin/smc2bed --log-file ${log_dir}/${arg_block}_out.log ${smc_dir}/${arg_block}_out.${iter}.smc.gz | bgzip > ${out_dir}/${arg_block}_out.${iter}.bed.gz 
#########################################################################################################

#########################################################################################################
# input arguments
args <- commandArgs(trailingOnly = TRUE)
arg_block   = args[1]
iter_start  = args[2]  
iter_end    = args[3]   
iter_jumps  = args[4]  # can take multiplications of 10 
smc_dir     = args[5]  # <-- SOURCE DIR CONTAINING SMC FILES (this is set in make-command_lines_for_beds.R)
log_dir     = args[6]  # <-- SOURCE DIR CONTAINING LOG FILES (this is set in make-command_lines_for_beds.R)
out_dir     = args[7]  # ./argBedFiles
log_file  = sprintf("%s/%s_out.log", log_dir , arg_block)
#########################################################################################################

## system(sprintf("cp %s %s", path_to_log_file, path_to_new_log))
#########################################################################################################

# create output directory if does not exists. Suppress warning (if dir already exists)
dir.create(out_dir,showWarnings=FALSE)

iterations = seq(as.numeric(iter_start),as.numeric(iter_end),as.numeric(iter_jumps))

#########################################################################################################
# main loop
for(iter in iterations) {
  smc_file     = sprintf("%s/%s_out.%s.smc.gz", smc_dir,arg_block, iter)
  out_bed_file = sprintf("%s/%s_out.%s.bed.gz", out_dir,arg_block, iter)
  cat("Creating bed file",out_bed_file,"\n")
  # make bed file
  system(sprintf("smc2bed --log-file %s %s | bgzip > %s", log_file ,smc_file, out_bed_file)) # <-- MAY NEED TO SET PATH TO SMC2BED
  # make .tbi index file
  system(sprintf("tabix %s", out_bed_file))

} # end of for(iter)
#########################################################################################################

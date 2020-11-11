#!/bin/bash
# script for running cross coalescence analysis for all peaks or control blocks in parallel

# directory where all log files are dumped
logDir=crossCoal-logs/
# in case script is run on FST peaks, should we also produce a stats file with average values per 20kb block
doFSTstatFile=true
inFile=$1       # infoTables/FST-table-with-center.txt    or   infoTables/control-blocks.txt
outFile=$2      # infoTables/cross_coals_fst_peaks.tsv    or   infoTables/cross_coals_control_blocks.tsv
regionType=$3   # fst or control (or any other string != "fst")
startRow=$4     # FST peak table starts at line 2 and control block table starts ato row 1
endRow=`cat $inFile | wc -l`  # can also manually truncate if you don't want to run on entire file
numInBatch=2   # number of analyses to be run in parallel. Set this based on your machine's capacity
scriptDir=`echo $0 | sed -r 's|^(.*)\/.*$|\1|'`

numJobs=0
tempOutFiles=""
rm -f $outFile
mkdir -p $logDir

for row in `seq $startRow $endRow`; do
   line=`head -n$row $inFile | tail -n1 | tr -d "\r" | tr -s "\t" " "`
   if [ "$regionType" == "fst" ]; then
      regionID=Contig`echo $line | cut -d" " -f1`
      scaffold=Contig`echo $line | cut -d" " -f2`
      peakStart=`     echo $line | cut -d" " -f4`
      peakEnd=`       echo $line | cut -d" " -f5`
      focusStart=`    echo $line | cut -d" " -f6`
      focusEnd=`      echo $line | cut -d" " -f7`
   else
      scaffold=`    echo $line | cut -d" " -f1`
      regionStart=` echo $line | cut -d" " -f2`
      regionEnd=`   echo $line | cut -d" " -f3`
      regionID=$scaffold.control_$((regionEnd/500000))
      peakStart=$((regionStart+150000))
      peakEnd=$((regionStart+350000))
      focusStart=$peakStart
      focusEnd=$peakEnd
   fi
   tempOut=`echo $outFile | sed -r 's/(.*)[.].*/\1/'`-$regionID.txt
   tempLog=$logDir/cross_coals_$regionID.log
   # echo "Rscript $scriptDir/cross_coal_per_region.R argTreeFiles/ $regionID $scaffold $peakStart $peakEnd $focusStart $focusEnd $tempOut"
   Rscript $scriptDir/cross_coal_per_region.R argTreeFiles/ $regionID $scaffold $peakStart $peakEnd $focusStart $focusEnd $tempOut &> $tempLog &
   tempOutFiles="$tempOutFiles $tempOut"
   ((numJobs++))
   # for Fst peaks, if doFSTstatFile flag is on, then create stats file
   if [ "$regionType" == "fst" ] && [ "$doFSTstatFile" == "true" ]; then
      tempLog=$logDir/cross_coals_${regionID}_stats.log
      Rscript $scriptDir/../smc_to_stats/window_cross_coal.R argTreeFiles/ $scaffold $((peakStart-200000)) $((peakEnd+200000)) 20000 argStats-windowed &> $tempLog &
      ((numJobs++)) 
   fi
   if [ $numJobs -eq $numInBatch ] || [ $row -eq $endRow ]; then
      wait
      if [ ! -e $outFile ]; then
         head -n1 -q $tempOutFiles | sort | uniq > $outFile
      fi
      tail -n1 -q $tempOutFiles >> $outFile
      if [ $? -ne 0 ]; then
         echo "Out files $tempOutFiles were not created. See log files in $logDir"
         exit 1
      fi
      rm -f $tempOutFiles
      numJobs=0
      tempOutFiles=""
   fi
done

type=$1		#0 or 1 - population index
dir=$2		#type: hard or soft or neutral
class=$3	#Directory: e.g. hardp2
windows=$4	#number of windows - 5
rep=$5		#number of replicates
windows=$(($windows-1)) #index of flanking window of interest
middle=$(($windows/2))	#index of middle window
path=$6 #path

for i in `seq 1 $rep`;
do
  for j in `seq 0 $windows`;
  do
    python local_LD.py $i $middle $j $type $dir $class $path
  done
done

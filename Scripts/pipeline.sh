type=$1		#Type of the simulation - hard or soft or neutral
rep=$2		#Number of replicates
dir=$3		#Location of the simulations - directory
id=$4		#id number
windows=$5	#Number of windows - 5
num_samples=$6	#Number of haploid samples - 48

for i in `seq 1 $rep`; do sh clean.sh $type$i.txt $dir; done #Clean the SLiM files by extracting the segsites, position, and genotypic data

##Global Statistics##
for i in `seq 1 $rep`; do python splitMsOutputIntoWindows.py $dir/_$type$i.txt $windows $dir/$type$i.txt; done #Split each simulation into 5 windows
for i in $dir/$type*msWin; do sh eval.sh $i; done;  #Computes summary statistics or features per window 
for i in `seq 1 $rep`; do sh global_combine.sh $i $type $dir $id; done; #Combines summary statistics across windows

##Local (population) Statistics##
cp parse_pop.py $dir/
cd $dir
for i in $type*msWin; do for j in `seq 0 1`; do python parse_pop.py $num_samples 1 $j $i; done; done #Parse population specific data
rm parse_pop.py
cd ..
for j in `seq 0 1`; do for i in $dir/$j'_'$type*msWin; do sh eval.sh $i; done; done #Computes summary statistics or features per window per population 
for j in `seq 0 1`; do for i in `seq 1 $rep`; do sh run.sh $i $j $type $dir $id; done; done #Combines summary statistics across windows per population
mkdir $dir/LD
for j in `seq 0 1`; #scan through population indices
do
  sh run_LD.sh $j $type $dir $windows $rep $path #Computes LD distribution per population
done

i=$1
j=$2
type=$3
dir=$4
id=$5
mkdir tmpStatDir$id
cp $dir/$j'_'$type$i.txt*stats tmpStatDir$id/;
python combineWinsIntoFeatureVec.py tmpStatDir$id/ > $dir/$j'_'$type$i.txt;
rm -f tmpStatDir$id/*;
rm -rf tmpStatDir$id

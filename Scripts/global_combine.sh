i=$1
type=$2
dir=$3
id=$4
mkdir tmpStatDir$id
cp $dir/$type$i.txt*stats tmpStatDir$id/;
python combineWinsIntoFeatureVecG.py tmpStatDir$id/ > $dir/g$type$i.txt;
rm -f tmpStatDir$id/*;
rm -rf tmpStatDir$id

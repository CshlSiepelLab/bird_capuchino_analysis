echo $1
x=$(echo "cat" $1)
eval $x | ./niceStats > $1.stats

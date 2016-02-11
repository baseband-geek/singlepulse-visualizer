#!/bin/bash


base=$1



infile="$base.singlepulse"
flagfile="$base.flag"
outfile=$base"_flagged.singlepulse"

while read line
do
    # load the flag_lo and flag_hi values from file
    flag_lo=`echo $line | awk '{print $1}'`
    flag_hi=`echo $line | awk '{print $2}'`
    timerange_list="$timerange_list !(\$3>$flag_lo && \$3<$flag_hi) &&"
done <$flagfile

echo $outfile

cat $infile | awk "{if ($timerange_list 1 ){print \$0;}}" > $outfile


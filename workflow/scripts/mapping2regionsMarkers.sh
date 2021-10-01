#!/usr/bin/env bash

set -e 

( echo "`date -R`: Sort the region file: $4..." &&
sort -u -k1,1 -k2,2n -k3,3n results/regions/$4 > $2/$3_$5_sorted.bed &&
echo "`date -R`: Success!" ||
{ echo "`date -R`: Process failed..."; exit 1; } ) > $9 2>&1 

( echo "`date -R`: Intersecting with $4..." &&
bedtools intersect \
-a $2/$3_$5_sorted.bed \
-b $1 \
-wa -c -F 0.5 > $2/$3_$5.txt &&
echo "`date -R`: Success!" ||
{ echo "`date -R`: Process failed..."; exit 1; } ) >> $9 2>&1

# no Windows
if [ $6 == "notCombine" ]; then

    ( echo "`date -R`: Passing combine..." &&
    cp $2/$3_$5.txt $2/$3_$5_combined.txt &&
    echo "`date -R`: Success!" ||
    { echo "`date -R`: Process failed..."; exit 1; } ) >> $9 2>&1

else

    if [ $6 == "-" ]; then 
        combOpt="" 
    else
        combOpt="$6"
    fi 

    (echo "`date -R`: Combine windows..." &&
    python3 workflow/scripts/combinewindows.py \
    -i $2/$3_$5.txt \
    $combOpt \
    -o $2/$3_$5_combined.txt &&
    echo "`date -R`: Success! Windows are combined." || 
    { echo "`date -R`: Process failed..."; exit 1; }  ) >> $9 2>&1

fi

moreinfo="$(echo $7 | sed 's/,/\t/g' | cut -f2 )"

( echo "`date -R`: Add info..." &&
python3 workflow/scripts/addColumns.py \
-i $2/$3_$5_combined.txt \
-o $2/$3_$5_combined_info.txt \
-c $moreinfo + $8 &&
echo "`date -R`: Success! Info is added." ||
{ echo "`date -R`: Process failed..."; exit 1; } ) >> $9 2>&1

(echo "`date -R`: Calculating RPKM values..." &&
python3 workflow/scripts/RPKM.py \
-i $2/$3_$5_combined_info.txt \
-chse 2 3 \
-c 7 \
-mr 0 \
-o $2/$3_$5_combined_rpkm.txt &&
echo "`date -R`: Success! RPKMs are calculated." || 
{ echo "`date -R`: Process failed..."; rm $2/$3_$5_combined_rpkm.txt; exit 1; }  ) >> $9 2>&1

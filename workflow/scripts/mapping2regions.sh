#!/usr/bin/env bash

set -e 

( echo "`date -R`: Sort the region file: $5..." &&
sort -u -k1,1 -k2,2n -k3,3n results/regions/$5 > $3/$4_$6_sorted.bed &&
echo "`date -R`: Success!" ||
{ echo "`date -R`: Process failed..."; exit 1; } ) > ${10} 2>&1 

( echo "`date -R`: Intersecting plus strand with $5..." &&
bedtools intersect \
-a $3/$4_$6_sorted.bed \
-b $1 \
-wa -c -F 0.5 > $3/$4_plus_$6.txt &&
echo "`date -R`: Success!" ||
{ echo "`date -R`: Process failed..."; exit 1; } ) >> ${10} 2>&1

( echo "`date -R`: Intersecting minus strand with $5..." &&
bedtools intersect \
-a $3/$4_$6_sorted.bed \
-b $2 \
-wa -c -F 0.5 > $3/$4_minus_$6.txt &&
echo "`date -R`: Success!" ||
{ echo "`date -R`: Process failed..."; exit 1; } ) >> ${10} 2>&1

# no Windows
if [ $7 == "notCombine" ]; then

    ( echo "`date -R`: Passing combine..." &&
    cp $3/$4_plus_$6.txt $3/$4_plus_$6_combined.txt &&
    cp $3/$4_minus_$6.txt $3/$4_minus_$6_combined.txt &&
    echo "`date -R`: Success!" ||
    { echo "`date -R`: Process failed..."; exit 1; } ) >> ${10} 2>&1

else

    if [ $7 == "-" ]; then 
        combOpt="" 
    else
        combOpt="$7"
    fi 

    (echo "`date -R`: Combine windows (plus strand)..." &&
    python3 workflow/scripts/combinewindows.py \
    -i $3/$4_plus_$6.txt \
    $combOpt \
    -o $3/$4_plus_$6_combined.txt &&
    echo "`date -R`: Success! Windows are combined." || 
    { echo "`date -R`: Process failed..."; exit 1; }  ) >> ${10} 2>&1

    (echo "`date -R`: Combine windows (minus strand)..." &&
    python3 workflow/scripts/combinewindows.py \
    -i $3/$4_minus_$6.txt \
    $combOpt \
    -o $3/$4_minus_$6_combined.txt &&
    echo "`date -R`: Success! Windows are combined." || 
    { echo "`date -R`: Process failed..."; exit 1; }  ) >> ${10} 2>&1

fi

moreinfo="$(echo $8 | sed 's/,/\t/g' )"

( echo "`date -R`: Add info (plus strand)..." &&
python3 workflow/scripts/addColumns.py \
-i $3/$4_plus_$6_combined.txt \
-o $3/$4_plus_$6_combined_info.txt \
-c $moreinfo + $9 &&
echo "`date -R`: Success! Info is added." ||
{ echo "`date -R`: Process failed..."; exit 1; } ) >> ${10} 2>&1

( echo "`date -R`: Add info (minus strand)..." &&
python3 workflow/scripts/addColumns.py \
-i $3/$4_minus_$6_combined.txt \
-o $3/$4_minus_$6_combined_info.txt \
-c $moreinfo - $9 &&
echo "`date -R`: Success! Info is added." ||
{ echo "`date -R`: Process failed..."; exit 1; } ) >> ${10} 2>&1

(echo "`date -R`: Calculating RPKM values (plus strand)..." &&
python3 workflow/scripts/RPKM.py \
-i $3/$4_plus_$6_combined_info.txt \
-chse 2 3 \
-c 7 \
-mr 0 \
-o $3/$4_plus_$6_combined_rpkm.txt &&
echo "`date -R`: Success! RPKMs are calculated." || 
{ echo "`date -R`: Process failed..."; rm $3/$4_plus_$6_combined_rpkm.txt; exit 1; }  ) >> ${10} 2>&1

(echo "`date -R`: Calculating RPKM values (minus strand)..." &&
python3 workflow/scripts/RPKM.py \
-i $3/$4_minus_$6_combined_info.txt \
-chse 2 3 \
-c 7 \
-mr 0 \
-o $3/$4_minus_$6_combined_rpkm.txt &&
echo "`date -R`: Success! RPKMs are calculated." || 
{ echo "`date -R`: Process failed..."; rm $3/$4_minus_$6_combined_rpkm.txt; exit 1; }  ) >> ${10} 2>&1

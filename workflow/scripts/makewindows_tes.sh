set -e

strands=("-" "+" ".")
names=(minus plus nostrand)

#

    awk -v a="$4" -v b="$5" -v n="$7" '{print $1"\t"int(($2+$3)/2-a/2-a*(b-1)/2)"\t"int(($2+$3)/2+a/2+a*(b-1)/2)"\t"n"\t"".""\t"$6}' $1 > $1$4$5temp1.bed
  
    python3 workflow/scripts/ExactMatchingIntervals.py -i $1$4$5temp1.bed -g $2 -o $1$4$5temp2.bed

    for i in $(seq 0 2); do

        awk -v strand="${strands[$i]}" '{if($6==strand){print}}' $1$4$5temp2.bed > $1$4$5temp2_${names[$i]}.bed 

        bedtools intersect -a $1$4$5temp2_${names[$i]}.bed -b $1$4$5temp2_${names[$i]}.bed -wa -wb | cut -f 1-6 | uniq -c | grep "^      1" | cut -c 9- > $1$4$5temp_unique_${names[$i]}.bed

        if [ "$6" = "reverse" ]; then # if reverse numbering of the windows needed

            bedtools makewindows -b $1$4$5temp_unique_${names[$i]}.bed -n $5 -i srcwinnum -reverse > $1$4$5temp_windowed_${names[$i]}.bed 

        else
            
            bedtools makewindows -b $1$4$5temp_unique_${names[$i]}.bed -n $5 -i srcwinnum > $1$4$5temp_windowed_${names[$i]}.bed 

        fi

        awk -v strand="${strands[$i]}" '{print $0"\t"".""\t"strand}' $1$4$5temp_windowed_${names[$i]}.bed > $1$4$5temp_windows_${names[$i]}.bed

    done
	
    cat $1$4$5temp_windows_${names[0]}.bed $1$4$5temp_windows_${names[1]}.bed $1$4$5temp_windows_${names[2]}.bed > $3

    

    [ -s $3 ] && echo "File not empty" || { echo "File empty"; rm $3; exit 1; }

    # Remove Unnecessary

        rm $1$4$5temp1.bed
        rm $1$4$5temp2.bed

        for name in ${names[@]}; do
        
            rm $1$4$5temp2_${name}.bed
            rm $1$4$5temp_unique_${name}.bed
            rm $1$4$5temp_windowed_${name}.bed
            rm $1$4$5temp_windows_${name}.bed

        done
    #

    


strands=("-" "+" ".")
names=(minus plus nostrand)

#

    awk -v a="$4" -v b="$5" '{print $1"\t"int(($2+$3)/2-a/2-a*(b-1)/2)"\t"int(($2+$3)/2+a/2+a*(b-1)/2)"\t"$4"\t"".""\t"$6}' $1 > temp1.bed
  
    python3 workflow/scripts/ExactMatchingIntervals.py -i temp1.bed -g $2 -o temp2.bed

    for i in $(seq 0 2); do

        awk -v strand="${strands[$i]}" '{if($6==strand){print}}' temp2.bed > temp2_${names[$i]}.bed 

        bedtools intersect -a temp2_${names[$i]}.bed -b temp2_${names[$i]}.bed -wa -wb | cut -f 1-6 | uniq -c | grep "^      1" | cut -c 9- > temp_unique_${names[$i]}.bed

        if [ "$6" = "reverse" ]; then # if reverse numbering of the windows needed

            bedtools makewindows -b temp_unique_${names[$i]}.bed -n $5 -i srcwinnum -reverse > temp_windowed_${names[$i]}.bed 

        else
            
            bedtools makewindows -b temp_unique_${names[$i]}.bed -n $5 -i srcwinnum > temp_windowed_${names[$i]}.bed 

        fi

        awk -v strand="${strands[$i]}" '{print $0"\t"".""\t"strand}' temp_windowed_${names[$i]}.bed > temp_windows_${names[$i]}.bed

    done
	
    cat temp_windows_${names[0]}.bed temp_windows_${names[1]}.bed temp_windows_${names[2]}.bed > $3
     
    # Remove Unnecessary

        rm temp1.bed
        rm temp2.bed

        for name in ${names[@]}; do
        
            rm temp2_${name}.bed
            rm temp_unique_${name}.bed
            rm temp_windowed_${name}.bed
            rm temp_windows_${name}.bed

        done
    #

    

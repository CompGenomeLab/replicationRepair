toolPath="/home/azgarian/Documents/myprojects/replicationRepair/gitignore/temp_codes"
dataPath="/home/azgarian/Desktop/masaustu"
genomePath="/home/azgarian/Documents/myprojects/replicationRepair/0_data/gitignore/Genome"
outputPath="/home/azgarian/Desktop"

intervalLen=100
windowNum=201
name="tss"

set -- "hg19_ucsc_genes_knownCanonical_minus_tss"

strands=("-" "+" ".")
names=(minus plus nostrand)

#

    awk -v a="$intervalLen" -v b="$windowNum" -v c="$name" '{print $1"\t"int(($2+$3)/2-a/2-a*(b-1)/2)"\t"int(($2+$3)/2+a/2+a*(b-1)/2)"\t"c"\t"".""\t"$6}' ${dataPath}/$1.bed > ${outputPath}/$1_organized2.bed
  
    ${toolPath}/ExactMatchingIntervals.py -i ${outputPath}/$1_organized2.bed -g ${genomePath}/genome_hg19.bed -o ${outputPath}/$1_organized.bed

    for i in $(seq 0 2); do

        awk -v strand="${strands[$i]}" '{if($6==strand){print}}' ${outputPath}/$1_organized.bed > ${outputPath}/$1_organized_${names[$i]}.bed 

        bedtools intersect -a ${outputPath}/$1_organized_${names[$i]}.bed -b ${outputPath}/$1_organized_${names[$i]}.bed -wa -wb | cut -f 1-6 | uniq -c | grep "^      1" | cut -c 9- > ${outputPath}/$1_unique_${names[$i]}.bed

        bedtools makewindows -b ${outputPath}/$1_unique_${names[$i]}.bed -n $windowNum -i srcwinnum -reverse > ${outputPath}/$1_windowed_${names[$i]}.bed 

        awk -v strand="${strands[$i]}" '{print $0"\t"".""\t"strand}' ${outputPath}/$1_windowed_${names[$i]}.bed > ${outputPath}/$1_windows_${names[$i]}.bed

    done
	
    cat ${outputPath}/$1_windows_${names[0]}.bed ${outputPath}/$1_windows_${names[1]}.bed ${outputPath}/$1_windows_${names[2]}.bed > ${outputPath}/$1_windows_${windowNum}_${intervalLen}.bed



     
     
    # Remove Unnecessary

        rm ${outputPath}/$1_organized2.bed
        rm ${outputPath}/$1_organized.bed

        for name in ${names[@]}; do
        
            rm ${outputPath}/$1_organized_${name}.bed
            rm ${outputPath}/$1_unique_${name}.bed
            rm ${outputPath}/$1_windowed_${name}.bed
            rm ${outputPath}/$1_windows_${name}.bed

        done
    #

    

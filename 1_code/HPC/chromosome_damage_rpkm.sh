#!/bin/bash

source ~/source_dir.sh

getname="$(echo $1 | sed 's/-/_/g')"
moreinfo="$(grep ${getname} ${mainPath}/project_repair_replication_all_samples.csv | sed 's/,/\t/g')"
if [ -f ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus_dipyrimidines.bed ]; then
    minus_line="$(grep -c '^' ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus_dipyrimidines.bed)"
    plus_line="$(grep -c '^' ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus_dipyrimidines.bed)"
    mappedReads=`echo "$minus_line + $plus_line" | bc` # In order to create mappedReads variable; "$1_cutadapt_sorted_plus/minus_dipyrimidines.bed" file should exist at the chosen directory or pre_analysis process must be done to create a new "$1_cutadapt_sorted_plus/minus_dipyrimidines.bed" file.
else 
    echo mappedReads does not exist!
fi
now=$(date +"[%Y.%m.%d]") # date

mkdir -p ${mainPath}/Damageseq/chromosomes
        
bedtools intersect -a ${genomePath}/genome_chr.bed -b ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus_dipyrimidines.bed -wa -c -F 0.5 >  ${mainPath}/Damageseq/chromosomes/$1_cutadapt_sorted_plus_chromosomes.txt

bedtools intersect -a ${genomePath}/genome_chr.bed -b ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus_dipyrimidines.bed -wa -c -F 0.5 >  ${mainPath}/Damageseq/chromosomes/$1_cutadapt_sorted_minus_chromosomes.txt

echo intersect is over!

${NGStoolkitPath}/addColumns.py -i ${mainPath}/Damageseq/chromosomes/$1_cutadapt_sorted_plus_chromosomes.txt -o ${mainPath}/Damageseq/chromosomes/$1_cutadapt_sorted_plus_chromosomes_full.txt -c "${moreinfo}" + "${mappedReads}"

${NGStoolkitPath}/addColumns.py -i ${mainPath}/Damageseq/chromosomes/$1_cutadapt_sorted_minus_chromosomes.txt -o ${mainPath}/Damageseq/chromosomes/$1_cutadapt_sorted_minus_chromosomes_full.txt -c "${moreinfo}" - "${mappedReads}"

echo info added!

python ${mainPath}/RPKM.py -i ${mainPath}/Damageseq/chromosomes/$1_cutadapt_sorted_plus_chromosomes_full.txt -c 7 -mr 0 -chse 2 3 -o ${mainPath}/Damageseq/chromosomes/$1_cutadapt_sorted_plus_chromosomes_rpkm.txt

python ${mainPath}/RPKM.py -i ${mainPath}/Damageseq/chromosomes/$1_cutadapt_sorted_minus_chromosomes_full.txt -c 7 -mr 0 -chse 2 3 -o ${mainPath}/Damageseq/chromosomes/$1_cutadapt_sorted_minus_chromosomes_rpkm.txt

echo RPKM calculated!

cat ${mainPath}/Damageseq/chromosomes/$1_cutadapt_sorted_plus_chromosomes_rpkm.txt >> ${mainPath}/Damageseq/chromosomes/${now}final_chromosome_damage.txt

cat ${mainPath}/Damageseq/chromosomes/$1_cutadapt_sorted_minus_chromosomes_rpkm.txt >> ${mainPath}/Damageseq/chromosomes/${now}final_chromosome_damage.txt

echo  data send to final!
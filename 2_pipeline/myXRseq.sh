#!/bin/bash

source ../1_code/source_dir.sh
source ${codePath}/source_key.sh
source ${codePath}/source_dataset.sh
source ${codePath}/functions_repairRep.sh

# Primary variables 

control="${mainPath}/3_output/XRseq/$1/control"
mkdir -p ${control} # directory for control files

preAnalysis="${mainPath}/3_output/XRseq/$1/pre_analysis"
mkdir -p ${preAnalysis} # directory for pre_analysis files

intersectCombine="${mainPath}/3_output/XRseq/$1/intersect_combine"
mkdir -p ${intersectCombine} # directory for intersected and combined files

final="${mainPath}/3_output/XRseq/$1/final"
mkdir -p ${final} # directory for final files

getname="$(echo $1 | sed 's/-/_/g')"
moreinfo="$(grep ${getname} ${mainPath}/0_data/samples.csv | sed 's/,/\t/g')"
if [ -f ${preAnalysis}/$1_cutadapt_sorted_chr.bed ]; then
    mappedReads="$(grep -c '^' ${preAnalysis}/$1_cutadapt_sorted_chr.bed)" # In order to create mappedReads variable; "$1_cutadapt_sorted_chr.bed" file should exist at the chosen directory or pre_analysis process must be done to create a new "$1_cutadapt_sorted_chr.bed" file.
else 
    echo mappedReads does not exist!
fi
now=$(date +"[%Y.%m.%d]") # date

# Pre-analysis

if ${Key_pre_analysis}; then

    check_raw_data $1 $rawdataPath  

    if ${Key_cutadapt}; then
    
        (cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG -o ${preAnalysis}/$1_cutadapt.fastq.gz ${rawdataPath}/$1${zip}) > ${control}/$1_control.txt 

        TotalFilteredReads="$(zcat ${preAnalysis}/$1_cutadapt.fastq.gz | echo $((`wc -l`/4)))"

        echo "Total Reads after cutadapt filtering: ${TotalFilteredReads}" >> ${control}/$1_control.txt 

        echo cutadapt is done!

    fi

    if ${Key_bowtie2}; then
    
        (bowtie2 -p 4 -x ${genomePath}/Bowtie2/genome -U ${preAnalysis}/$1_cutadapt.fastq.gz -S ${preAnalysis}/$1_cutadapt.sam) 2>> ${control}/$1_control.txt

        samtools view -Sb -o ${preAnalysis}/$1_cutadapt.bam ${preAnalysis}/$1_cutadapt.sam # samtools: sam to bam

        samtools view -q 20 -b ${preAnalysis}/$1_cutadapt.bam | bedtools bamtobed > ${preAnalysis}/$1_cutadapt.bed # bedtools: bam to bed

        rm ${preAnalysis}/$1_cutadapt.sam
        rm ${preAnalysis}/$1_cutadapt.fastq.gz

        echo bowtie2 mapping is done!

    fi

    if ${Key_sort_count}; then

        sort -u -k1,1 -k2,2n -k3,3n ${preAnalysis}/$1_cutadapt.bed > ${preAnalysis}/$1_cutadapt_sorted.bed # sort

        awk '{print $3-$2}' ${preAnalysis}/$1_cutadapt_sorted.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}' > ${control}/$1_cutadapt_length_distribution.txt # length distribution

        awk '{ if ($3-$2 == 26) { print } }' ${preAnalysis}/$1_cutadapt_sorted.bed > ${preAnalysis}/$1_cutadapt_sorted_26.bed # get 26 nucleotides long reads

        bedtools getfasta -fi ${genomePath}/genome.fa -bed ${preAnalysis}/$1_cutadapt_sorted_26.bed -fo ${preAnalysis}/$1_cutadapt_sorted_26.fa -s # bedtools: to FASTA format

        ${codePath}/fa2kmerAbundanceTable.py -i ${preAnalysis}/$1_cutadapt_sorted_26.fa -k 2 -o ${control}/$1_cutadapt_sorted_26_dinucleotideTable.txt # dinucleotide content of sequences at length 26

    fi

    if ${Key_sep_plus_minus}; then
    
        grep "chr" ${preAnalysis}/$1_cutadapt_sorted.bed | grep -v -e "chrY" -e "chrM" > ${preAnalysis}/$1_cutadapt_sorted_chr.bed

        mappedReads="$(grep -c '^' ${preAnalysis}/$1_cutadapt_sorted_chr.bed)"

        echo "Total Mapped Reads: ${mappedReads}" >> ${control}/$1_control.txt # count

        awk '{if($6=="+"){print}}' ${preAnalysis}/$1_cutadapt_sorted_chr.bed > ${preAnalysis}/$1_cutadapt_sorted_plus.bed # separating plus strand

        awk '{if($6=="-"){print}}' ${preAnalysis}/$1_cutadapt_sorted_chr.bed > ${preAnalysis}/$1_cutadapt_sorted_minus.bed # separating minus strand

        echo bed files are ready!

    fi

    if ${Key_bedgraph_BigWig}; then
    
        bedtools genomecov -i ${preAnalysis}/$1_cutadapt_sorted_plus.bed -g ${genomePath}/genome.fa.fai -bg -scale $(echo ${mappedReads} | awk '{print 1000000/$1}') > ${preAnalysis}/$1_cutadapt_sorted_plus.bdg # bedtools: to generate bedgraph files

        bedtools genomecov -i ${preAnalysis}/$1_cutadapt_sorted_minus.bed -g ${genomePath}/genome.fa.fai -bg -scale $(echo ${mappedReads} | awk '{print -1000000/$1}') > ${preAnalysis}/$1_cutadapt_sorted_minus.bdg # bedtools: to generate bedgraph files

        bedGraphToBigWig ${preAnalysis}/$1_cutadapt_sorted_plus.bdg ${genomePath}/genome.fa.fai ${preAnalysis}/$1_cutadapt_sorted_plus.bw # bedgraph to BigWig

        bedGraphToBigWig ${preAnalysis}/$1_cutadapt_sorted_minus.bdg ${genomePath}/genome.fa.fai ${preAnalysis}/$1_cutadapt_sorted_minus.bw # bedgraph to BigWig

        echo BigWig files are ready!

    fi

    if ${Key_TS_NTS}; then

        bedtools intersect -sorted -a ${genomePath}/hg19_ucsc_genes_knownCanonical_stranded.bed -b ${preAnalysis}/$1_cutadapt_sorted_chr.bed -wa -c -S -F 0.5 > ${preAnalysis}/$1_cutadapt_sorted_TScount.txt

        bedtools intersect -sorted -a ${genomePath}/hg19_ucsc_genes_knownCanonical_stranded.bed -b ${preAnalysis}/$1_cutadapt_sorted_chr.bed -wa -c -s -F 0.5 > ${preAnalysis}/$1_cutadapt_sorted_NTScount.txt

        paste ${preAnalysis}/$1_cutadapt_sorted_TScount.txt ${preAnalysis}/$1_cutadapt_sorted_NTScount.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$14}' > ${preAnalysis}/$1_cutadapt_sorted_TSoverNTScount.txt

        echo TS/NTS counts are calculated!
    
    fi

fi


# Analysis

if ${Key_downstream_analysis}; then

    for ((i=0;i<${#dataset[@]};i++)); do

        echo data: ${dataset[i]} long name: ${data_name[i]} combine options: ${combine_options[i]}

        if ${Key_alignment}; then

            if ${Key_intergenic}; then

                bedtools intersect -a ${preAnalysis}/$1_cutadapt_sorted_plus.bed -b ${genomePath}/hg19_ucsc_genes_knownCanonical_stranded.bed -v -f 0.5 >  ${intersectCombine}/$1_cutadapt_sorted_plus_filtered.bed

                bedtools intersect -a ${preAnalysis}/$1_cutadapt_sorted_minus.bed -b ${genomePath}/hg19_ucsc_genes_knownCanonical_stranded.bed -v -f 0.5 >  ${intersectCombine}/$1_cutadapt_sorted_minus_filtered.bed

                plus_line="$(grep -c '^' ${intersectCombine}/$1_cutadapt_sorted_plus_filtered.bed)"
                minus_line="$(grep -c '^' ${intersectCombine}/$1_cutadapt_sorted_minus_filtered.bed)"
                mappedReads=`echo "$minus_line + $plus_line" | bc`

                bedtools intersect -a ${mainPath}/data/${data_name[i]} -b ${intersectCombine}/$1_cutadapt_sorted_plus_filtered.bed -wa -c -F 0.5 >  ${intersectCombine}/$1_cutadapt_sorted_plus_${dataset[i]}.txt

                bedtools intersect -a ${mainPath}/data/${data_name[i]} -b ${intersectCombine}/$1_cutadapt_sorted_minus_filtered.bed -wa -c -F 0.5 >  ${intersectCombine}/$1_cutadapt_sorted_minus_${dataset[i]}.txt

            else

                bedtools intersect -a ${mainPath}/data/${data_name[i]} -b ${preAnalysis}/$1_cutadapt_sorted_plus.bed -wa -c -F 0.5 >  ${intersectCombine}/$1_cutadapt_sorted_plus_${dataset[i]}.txt

                bedtools intersect -a ${mainPath}/data/${data_name[i]} -b ${preAnalysis}/$1_cutadapt_sorted_minus.bed -wa -c -F 0.5 >  ${intersectCombine}/$1_cutadapt_sorted_minus_${dataset[i]}.txt

            fi

            echo ${dataset[i]} alignment done!
    
        fi

        if [[ ${dataset[i]} == *"windows"* ]]; then
            
            if ${Key_combineWindows}; then

                ${mainPath}/combinewindows.py -i ${intersectCombine}/$1_cutadapt_sorted_plus_${dataset[i]}.txt ${combine_options[i]} -o ${intersectCombine}/$1_cutadapt_sorted_plus_${dataset[i]}_combined.txt

                ${mainPath}/combinewindows.py -i ${intersectCombine}/$1_cutadapt_sorted_minus_${dataset[i]}.txt ${combine_options[i]} -o ${intersectCombine}/$1_cutadapt_sorted_minus_${dataset[i]}_combined.txt

                echo ${dataset[i]} combined!

            fi
        
            if ${Key_moreInfo}; then

                ${codePath}/addColumns.py -i ${intersectCombine}/$1_cutadapt_sorted_plus_${dataset[i]}_combined.txt -o ${final}/$1_cutadapt_sorted_plus_${dataset[i]}_full.txt -c "${moreinfo}" + "${mappedReads}"

                ${codePath}/addColumns.py -i ${intersectCombine}/$1_cutadapt_sorted_minus_${dataset[i]}_combined.txt -o ${final}/$1_cutadapt_sorted_minus_${dataset[i]}_full.txt -c "${moreinfo}" - "${mappedReads}"

                echo ${dataset[i]} info added!

            fi

        else 

            if ${Key_moreInfo}; then

                ${codePath}/addColumns.py -i ${intersectCombine}/$1_cutadapt_sorted_plus_${dataset[i]}.txt -o ${final}/$1_cutadapt_sorted_plus_${dataset[i]}_full.txt -c "${moreinfo}" + "${mappedReads}"

                ${codePath}/addColumns.py -i ${intersectCombine}/$1_cutadapt_sorted_minus_${dataset[i]}.txt -o ${final}/$1_cutadapt_sorted_minus_${dataset[i]}_full.txt -c "${moreinfo}" - "${mappedReads}"

                echo ${dataset[i]} info added!

            fi

        fi

        if ${Key_rpkm}; then
        
            python ${mainPath}/RPKM.py -i ${final}/$1_cutadapt_sorted_plus_${dataset[i]}_full.txt -chse 2 3 -c 7 -mr 0 -o ${final}/$1_cutadapt_sorted_plus_${dataset[i]}_rpkm.txt

            python ${mainPath}/RPKM.py -i ${final}/$1_cutadapt_sorted_minus_${dataset[i]}_full.txt -chse 2 3 -c 7 -mr 0 -o ${final}/$1_cutadapt_sorted_minus_${dataset[i]}_rpkm.txt

            echo ${dataset[i]} rpkm calculated!

        fi

        if ${Key_file}; then

            cat ${final}/$1_cutadapt_sorted_plus_${dataset[i]}_rpkm.txt >> ${mainPath}/${now}final_report_${dataset[i]}.txt

            cat ${final}/$1_cutadapt_sorted_minus_${dataset[i]}_rpkm.txt >> ${mainPath}/${now}final_report_${dataset[i]}.txt 

            echo ${dataset[i]} final report created!

        fi

    done
    
fi

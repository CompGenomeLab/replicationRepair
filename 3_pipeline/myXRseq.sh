#!/bin/bash

source ~/repairReplication/source_dir.sh
source ${mainPath}/source_key.sh
source ${mainPath}/source_dataset.sh
source ${mainPath}/functions_repairRep.sh

# Primary variables 

    getname="$(echo $1 | sed 's/-/_/g')"
    moreinfo="$(grep ${getname} ${mainPath}/project_repair_replication_all_samples.csv | sed 's/,/\t/g')"
    if [ -f ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed ]; then
        mappedReads="$(grep -c '^' ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed)" # In order to create mappedReads variable; "$1_cutadapt_sorted_chr.bed" file should exist at the chosen directory or pre_analysis process must be done to create a new "$1_cutadapt_sorted_chr.bed" file.
    else 
        echo mappedReads does not exist!
    fi
    now=$(date +"[%Y.%m.%d]") # date

    mkdir -p ${mainPath}/XRseq/$1/control # directory for control files
    mkdir -p ${mainPath}/XRseq/$1/pre_analysis # directory for pre_analysis files
    mkdir -p ${mainPath}/XRseq/$1/intersect_combine # directory for intersected and combined files
    mkdir -p ${mainPath}/XRseq/$1/final # directory for final files

# Pre-analysis

    if ${Key_pre_analysis}; then

        check_raw_data $1 $rawdataPath  

        TotalFastqReads="$(zcat ${rawdataPath}/$1${zip} | echo $((`wc -l`/4)))"
	
	    echo "Total Fastq Reads: ${TotalFastqReads}" > ${mainPath}/XRseq/$1/control/$1_control.txt # Control: add total fastq reads

        if ${Key_cutadapt}; then
        
            (cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG -o ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt.fastq.gz ${rawdataPath}/$1${zip}) >> ${mainPath}/XRseq/$1/control/$1_control.txt 

	        TotalFilteredReads="$(zcat ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt.fastq.gz | echo $((`wc -l`/4)))"

	        echo "Total Reads after cutadapt filtering: ${TotalFilteredReads}" >> ${mainPath}/XRseq/$1/control/$1_control.txt 

            echo cutadapt is done!

        fi

        if ${Key_bowtie2}; then
        
            (bowtie2 -p 4 -x ${genomePath}/Bowtie2/genome -U ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt.fastq.gz -S ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt.sam) 2>> ${mainPath}/XRseq/$1/control/$1_control.txt

	        samtools view -Sb -o ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt.bam ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt.sam # samtools: sam to bam

	        samtools view -q 20 -b ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt.bam | bedtools bamtobed > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt.bed # bedtools: bam to bed

            rm ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt.sam
            rm ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt.fastq.gz

            echo bowtie2 mapping is done!

        fi

        if ${Key_sort_count}; then

	        sort -u -k1,1 -k2,2n -k3,3n ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt.bed > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted.bed # sort

	        awk '{print $3-$2}' ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}' > ${mainPath}/XRseq/$1/control/$1_cutadapt_length_distribution.txt # length distribution

	        awk '{ if ($3-$2 == 26) { print } }' ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted.bed > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_26.bed # get 26 nucleotides long reads

	        bedtools getfasta -fi ${genomePath}/genome.fa -bed ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_26.bed -fo ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_26.fa # bedtools: to FASTA format

	        ${NGStoolkitPath}/fa2kmerAbundanceTable.py -i ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_26.fa -k 2 -o ${mainPath}/XRseq/$1/control/$1_cutadapt_sorted_26_dinucleotideTable.txt # dinucleotide content of sequences at length 26

        fi

        if ${Key_sep_plus_minus}; then
        
	        grep "chr" ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted.bed | grep -v -e "chrY" -e "chrM" > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed

            mappedReads="$(grep -c '^' ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed)"

            echo "Total Mapped Reads: ${mappedReads}" >> ${mainPath}/XRseq/$1/control/$1_control.txt # count

	        awk '{if($6=="+"){print}}' ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_plus.bed # separating plus strand

	        awk '{if($6=="-"){print}}' ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_minus.bed # separating minus strand

            echo bed files are ready!

        fi

        if ${Key_bedgraph_BigWig}; then
        
	        bedtools genomecov -i ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_plus.bed -g ${genomePath}/genome.fa.fai -bg -scale $(echo ${mappedReads} | awk '{print 1000000/$1}') > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_plus.bdg # bedtools: to generate bedgraph files

	        bedtools genomecov -i ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_minus.bed -g ${genomePath}/genome.fa.fai -bg -scale $(echo ${mappedReads} | awk '{print -1000000/$1}') > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_minus.bdg # bedtools: to generate bedgraph files

	        bedGraphToBigWig ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_plus.bdg ${genomePath}/genome.fa.fai ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_plus.bw # bedgraph to BigWig

	        bedGraphToBigWig ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_minus.bdg ${genomePath}/genome.fa.fai ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_minus.bw # bedgraph to BigWig

            echo BigWig files are ready!

        fi

        if ${Key_TS_NTS}; then

            bedtools intersect -sorted -a ${genomePath}/ensembl_genes.bed -b ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed -wa -c -S -F 0.5 > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TScount.txt

            bedtools intersect -sorted -a ${genomePath}/ensembl_genes.bed -b ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed -wa -c -s -F 0.5 > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_NTScount.txt

            paste ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TScount.txt ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_NTScount.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$14}' > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TSoverNTScount.txt

            echo TS/NTS counts are calculated!
        
        fi
    
    fi


# Analysis

    if ${Key_downstream_analysis}; then

        for ((i=0;i<${#dataset[@]};i++)); do

            echo data: ${dataset[i]} long name: ${data_name[i]} combine options: ${combine_options[i]}

            if ${Key_alignment}; then


                bedtools intersect -a ${mainPath}/data/${data_name[i]} -b ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_plus.bed -wa -c -F 0.5 >  ${mainPath}/XRseq/$1/intersect_combine/$1_cutadapt_sorted_plus_${dataset[i]}.txt

                bedtools intersect -a ${mainPath}/data/${data_name[i]} -b ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_minus.bed -wa -c -F 0.5 >  ${mainPath}/XRseq/$1/intersect_combine/$1_cutadapt_sorted_minus_${dataset[i]}.txt

                echo ${dataset[i]} alignment done!
        
            fi

            if [[ ${dataset[i]} == *"windows"* ]]; then
                
                if ${Key_combineWindows}; then

		            ${mainPath}/combinewindows.py -i ${mainPath}/XRseq/$1/intersect_combine/$1_cutadapt_sorted_plus_${dataset[i]}.txt ${combine_options[i]} -o ${mainPath}/XRseq/$1/intersect_combine/$1_cutadapt_sorted_plus_${dataset[i]}_combined.txt

		            ${mainPath}/combinewindows.py -i ${mainPath}/XRseq/$1/intersect_combine/$1_cutadapt_sorted_minus_${dataset[i]}.txt ${combine_options[i]} -o ${mainPath}/XRseq/$1/intersect_combine/$1_cutadapt_sorted_minus_${dataset[i]}_combined.txt

                    echo ${dataset[i]} combined!

                fi
            
                if ${Key_moreInfo}; then

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/XRseq/$1/intersect_combine/$1_cutadapt_sorted_plus_${dataset[i]}_combined.txt -o ${mainPath}/XRseq/$1/final/$1_cutadapt_sorted_plus_${dataset[i]}_full.txt -c "${moreinfo}" + "${mappedReads}"

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/XRseq/$1/intersect_combine/$1_cutadapt_sorted_minus_${dataset[i]}_combined.txt -o ${mainPath}/XRseq/$1/final/$1_cutadapt_sorted_minus_${dataset[i]}_full.txt -c "${moreinfo}" - "${mappedReads}"

                    echo ${dataset[i]} info added!

                fi

            else 

                if ${Key_moreInfo}; then

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/XRseq/$1/intersect_combine/$1_cutadapt_sorted_plus_${dataset[i]}.txt -o ${mainPath}/XRseq/$1/final/$1_cutadapt_sorted_plus_${dataset[i]}_full.txt -c "${moreinfo}" + "${mappedReads}"

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/XRseq/$1/intersect_combine/$1_cutadapt_sorted_minus_${dataset[i]}.txt -o ${mainPath}/XRseq/$1/final/$1_cutadapt_sorted_minus_${dataset[i]}_full.txt -c "${moreinfo}" - "${mappedReads}"

                    echo ${dataset[i]} info added!

                fi

            fi

            if ${Key_rpkm}; then
            
                python ${mainPath}/RPKM.py -i ${mainPath}/XRseq/$1/final/$1_cutadapt_sorted_plus_${dataset[i]}_full.txt -chse 2 3 -c 7 -mr 0 -o ${mainPath}/XRseq/$1/final/$1_cutadapt_sorted_plus_${dataset[i]}_rpkm.txt

                python ${mainPath}/RPKM.py -i ${mainPath}/XRseq/$1/final/$1_cutadapt_sorted_minus_${dataset[i]}_full.txt -chse 2 3 -c 7 -mr 0 -o ${mainPath}/XRseq/$1/final/$1_cutadapt_sorted_minus_${dataset[i]}_rpkm.txt

                echo ${dataset[i]} rpkm calculated!

            fi

            if ${Key_file}; then

                cat ${mainPath}/XRseq/$1/final/$1_cutadapt_sorted_plus_${dataset[i]}_rpkm.txt >> ${mainPath}/${now}final_report_${dataset[i]}.txt

                cat ${mainPath}/XRseq/$1/final/$1_cutadapt_sorted_minus_${dataset[i]}_rpkm.txt >> ${mainPath}/${now}final_report_${dataset[i]}.txt 

                echo ${dataset[i]} final report created!

            fi

        done
        
    fi

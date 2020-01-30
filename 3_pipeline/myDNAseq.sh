#!/bin/bash

source ~/repairReplication/source_dir.sh
source ${mainPath}/source_key.sh
source ${mainPath}/source_dataset.sh
source ${mainPath}/functions_repairRep.sh

# Primary variables 

    getname="$(echo $1 | sed 's/-/_/g')"
    moreinfo="$(grep ${getname} ${mainPath}/project_repair_replication_all_samples.csv | sed 's/,/\t/g')"
    if [ -f ${mainPath}/Input/$1/pre_analysis/$1_chr.bed ]; then
    	mappedReads="$(grep -c '^' ${mainPath}/Input/$1/pre_analysis/$1_chr.bed)" # In order to create mappedReads variable; "$1_chr.bed" file should exist at the chosen directory or pre_analysis process must be done to create a new "$1_chr.bed" file.

    else 
        echo mappedReads does not exist!
    fi

    now=$(date +"[%Y.%m.%d]") # date

    mkdir -p ${mainPath}/Input/$1/control # directory for control files
    mkdir -p ${mainPath}/Input/$1/pre_analysis # directory for pre_analysis files
    mkdir -p ${mainPath}/Input/$1/intersect_combine # directory for intersected and combined files
    mkdir -p ${mainPath}/Input/$1/final # directory for final files


# Pre-analysis

    if ${Key_pre_analysis}; then

        check_raw_data $1 $rawdataPath    
        
        if [ $layout == "p" ]; then

            TotalFastqReads_R1="$(grep -c '@' ${rawdataPath}/$1${zip})"

            TotalFastqReads_R2="$(grep -c '@' ${rawdataPath}/$1${zip})"

            TotalFastqReads="$(($TotalFastqReads_R1+$TotalFastqReads_R2))"

	        echo "Total Fastq Reads: ${TotalFastqReads}" > ${mainPath}/Input/$1/control/$1_control.txt # Control: add total fastq reads

            if ${Key_bowtie2}; then
            
	            (bowtie2 -p 4 -X 1000 -x ${genomePath}/Bowtie2/genome -1 ${rawdataPath}/$1${zip} -2 ${rawdataPath}/$1${zip} -S ${mainPath}/Input/$1/pre_analysis/$1.sam) 2>> ${mainPath}/Input/$1/control/$1_control.txt

	            samtools view -Sb -o ${mainPath}/Input/$1/pre_analysis/$1_cutadapt.bam ${mainPath}/Input/$1/pre_analysis/$1_cutadapt.sam # samtools: sam to bam

                samtools view -q 20 -bf 0x2 ${mainPath}/Input/$1/pre_analysis/$1_cutadapt.bam | sort -n | bedtools bamtobed -bedpe -mate1 > ${mainPath}/Input/$1/pre_analysis/$1_cutadapt.bedpe

                awk '{if($9=="+"){print $1"\t"$2"\t"$6"\t"$7"\t"$8"\t"$9}}' ${mainPath}/Input/$1/pre_analysis/$1_cutadapt.bedpe > ${mainPath}/Input/$1/pre_analysis/$1_cutadapt_plus.bed

                awk '{if($9=="-"){print $1"\t"$5"\t"$3"\t"$7"\t"$8"\t"$9}}' ${mainPath}/Input/$1/pre_analysis/$1_cutadapt.bedpe > ${mainPath}/Input/$1/pre_analysis/$1_cutadapt_minus.bed
                
                cat ${mainPath}/Input/$1/pre_analysis/$1_cutadapt_plus.bed ${mainPath}/Input/$1/pre_analysis/$1_cutadapt_minus.bed > ${mainPath}/Input/$1/pre_analysis/$1_cutadapt.bed
                
                rm ${mainPath}/Input/$1/pre_analysis/$1.sam

                echo bowtie2 mapping is done!

            fi
        
        elif [ $layout == "s" ]; then

            TotalFastqReads="$(grep -c '@' ${rawdataPath}/$1${zip})"
        
            echo "Total Fastq Reads: ${TotalFastqReads}" > ${mainPath}/Input/$1/control/$1_control.txt # Control: add total fastq reads

            if ${Key_bowtie2}; then
            
                (bowtie2 -p 4 -x ${genomePath}/Bowtie2/genome -U ${rawdataPath}/$1${zip} -S ${mainPath}/Input/$1/pre_analysis/$1.sam) 2>> ${mainPath}/Input/$1/control/$1_control.txt

                samtools view -Sb -o ${mainPath}/Input/$1/pre_analysis/$1.bam ${mainPath}/Input/$1/pre_analysis/$1.sam # samtools: sam to bam

                samtools view -q 20 -b ${mainPath}/Input/$1/pre_analysis/$1.bam | bedtools bamtobed  > ${mainPath}/Input/$1/pre_analysis/$1.bed # bedtools: bam to bed

                rm ${mainPath}/Input/$1/pre_analysis/$1.sam

                echo bowtie2 mapping is done!

            fi

        fi

        if ${Key_sort_count}; then

	        sort -u -k1,1 -k2,2n -k3,3n ${mainPath}/Input/$1/pre_analysis/$1.bed > ${mainPath}/Input/$1/pre_analysis/$1_sorted.bed # sort

	        awk '{print $3-$2}' ${mainPath}/Input/$1/pre_analysis/$1_sorted.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}' > ${mainPath}/Input/$1/control/$1_length_distribution.txt # length distribution

        fi

        if ${Key_sep_plus_minus}; then
        
	        grep "chr" ${mainPath}/Input/$1/pre_analysis/$1_sorted.bed | grep -v -e "chrY" -e "chrM" > ${mainPath}/Input/$1/pre_analysis/$1_chr.bed

            mappedReads="$(grep -c '^' ${mainPath}/Input/$1/pre_analysis/$1_chr.bed)"

            echo "Total Mapped Reads: ${mappedReads}" >> ${mainPath}/Input/$1/control/$1_control.txt # count

	        awk '{if($6=="+"){print}}' ${mainPath}/Input/$1/pre_analysis/$1_chr.bed > ${mainPath}/Input/$1/pre_analysis/$1_plus.bed # separating plus strand

	        awk '{if($6=="-"){print}}' ${mainPath}/Input/$1/pre_analysis/$1_chr.bed > ${mainPath}/Input/$1/pre_analysis/$1_minus.bed # separating minus strand

            echo bed files are ready!

        fi

          
    fi


# Analysis

    if ${Key_downstream_analysis}; then

        for ((i=0;i<${#dataset[@]};i++)); do

            echo data: ${dataset[i]} long name: ${data_name[i]} combine options: ${combine_options[i]}

            if ${Key_alignment}; then

                bedtools intersect -a ${mainPath}/data/${data_name[i]} -b ${mainPath}/Input/$1/pre_analysis/$1_plus.bed -wa -c -F 0.5 >  ${mainPath}/Input/$1/intersect_combine/$1_plus_${dataset[i]}.txt

                bedtools intersect -a ${mainPath}/data/${data_name[i]} -b ${mainPath}/Input/$1/pre_analysis/$1_minus.bed -wa -c -F 0.5 >  ${mainPath}/Input/$1/intersect_combine/$1_minus_${dataset[i]}.txt

                echo ${dataset[i]} alignment done!
        
            fi

            if [[ ${dataset[i]} == *"windows"* ]]; then
                
                if ${Key_combineWindows}; then

		            ${mainPath}/combinewindows.py -i ${mainPath}/Input/$1/intersect_combine/$1_plus_${dataset[i]}.txt ${combine_options[i]} -o ${mainPath}/Input/$1/intersect_combine/$1_plus_${dataset[i]}_combined.txt

		            ${mainPath}/combinewindows.py -i ${mainPath}/Input/$1/intersect_combine/$1_minus_${dataset[i]}.txt ${combine_options[i]} -o ${mainPath}/Input/$1/intersect_combine/$1_minus_${dataset[i]}_combined.txt

                    echo ${dataset[i]} combined!

                fi
            
                if ${Key_moreInfo}; then

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/Input/$1/intersect_combine/$1_plus_${dataset[i]}_combined.txt -o ${mainPath}/Input/$1/final/$1_plus_${dataset[i]}_full.txt -c "${moreinfo}" + "${mappedReads}"

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/Input/$1/intersect_combine/$1_minus_${dataset[i]}_combined.txt -o ${mainPath}/Input/$1/final/$1_minus_${dataset[i]}_full.txt -c "${moreinfo}" - "${mappedReads}"

                    echo ${dataset[i]} info added!

                fi

            else 

                if ${Key_moreInfo}; then

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/Input/$1/intersect_combine/$1_plus_${dataset[i]}.txt -o ${mainPath}/Input/$1/final/$1_plus_${dataset[i]}_full.txt -c "${moreinfo}" + "${mappedReads}"

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/Input/$1/intersect_combine/$1_minus_${dataset[i]}.txt -o ${mainPath}/Input/$1/final/$1_minus_${dataset[i]}_full.txt -c "${moreinfo}" - "${mappedReads}"

                    echo ${dataset[i]} info added!

                fi

            fi

            if ${Key_rpkm}; then
            
                python ${mainPath}/RPKM.py -i ${mainPath}/Input/$1/final/$1_plus_${dataset[i]}_full.txt -chse 2 3 -c 7 -mr 0 -o ${mainPath}/Input/$1/final/$1_plus_${dataset[i]}_rpkm.txt

                python ${mainPath}/RPKM.py -i ${mainPath}/Input/$1/final/$1_minus_${dataset[i]}_full.txt -chse 2 3 -c 7 -mr 0 -o ${mainPath}/Input/$1/final/$1_minus_${dataset[i]}_rpkm.txt

                echo ${dataset[i]} rpkm calculated!

            fi

            if ${Key_file}; then

                cat ${mainPath}/Input/$1/final/$1_plus_${dataset[i]}_rpkm.txt >> ${mainPath}/${now}final_report_${dataset[i]}.txt

                cat ${mainPath}/Input/$1/final/$1_minus_${dataset[i]}_rpkm.txt >> ${mainPath}/${now}final_report_${dataset[i]}.txt 

                echo ${dataset[i]} final report created!

            fi

        done
        
    fi



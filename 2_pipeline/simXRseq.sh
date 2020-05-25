#!/bin/bash

source ~/repairReplication/source_dir.sh
source ${mainPath}/source_key.sh
source ${mainPath}/source_dataset.sh
source ${mainPath}/functions_repairRep.sh

# Primary variables 

    getname="$(echo $1 | sed 's/-/_/g')"
    moreinfo="$(grep ${getname} ${mainPath}/project_repair_replication_all_samples.csv | sed 's/,/\t/g')"
    if [ -f ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_chr.bed ]; then
        mappedReads="$(grep -c '^' ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_chr.bed)" # In order to create mappedReads variable; "$1_sorted_chr.bed" file should exist at the chosen directory or pre_analysis process must be done to create a new "$1_sorted_chr.bed" file.
    else 
        echo mappedReads does not exist!
    fi
    now=$(date +"[%Y.%m.%d]") # date

    mkdir -p ${mainPath}/simulation/XRseq/$1/control # directory for control files
    mkdir -p ${mainPath}/simulation/XRseq/$1/pre_analysis # directory for pre_analysis files
    mkdir -p ${mainPath}/simulation/XRseq/$1/intersect_combine # directory for intersected and combined files
    mkdir -p ${mainPath}/simulation/XRseq/$1/final # directory for final files

# Pre-analysis

    if ${Key_pre_analysis}; then

        check_raw_data $1_sim $simRawPath 

        TotalFastqReads="$(grep -c '@' ${$simRawPath}/$1_sim${zip})"
	
	    echo "Total Fastq Reads: ${TotalFastqReads}" > ${mainPath}/simulation/XRseq/$1/control/$1_control.txt # Control: add total fastq reads

        if ${Key_bowtie2}; then
        
            (bowtie2 -p 4 -x ${genomePath}/Bowtie2/genome -U ${simRawPath}/$1_sim${zip} -S ${mainPath}/simulation/XRseq/$1/pre_analysis/$1.sam) 2>> ${mainPath}/simulation/XRseq/$1/control/$1_control.txt

	        samtools view -q 20 -Sb -o ${mainPath}/simulation/XRseq/$1/pre_analysis/$1.bam ${mainPath}/simulation/XRseq/$1/pre_analysis/$1.sam # samtools: sam to bam

	        bedtools bamtobed -i ${mainPath}/simulation/XRseq/$1/pre_analysis/$1.bam > ${mainPath}/simulation/XRseq/$1/pre_analysis/$1.bed # bedtools: bam to bed

            rm ${mainPath}/simulation/XRseq/$1/pre_analysis/$1.sam
            rm ${mainPath}/simulation/XRseq/$1/pre_analysis/$1.fastq.gz

            echo bowtie2 mapping is done!

        fi

        if ${Key_sort_count}; then

	        sort -u -k1,1 -k2,2n -k3,3n ${mainPath}/simulation/XRseq/$1/pre_analysis/$1.bed > ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted.bed # sort

	        bedtools getfasta -fi ${genomePath}/genome.fa -bed ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted.bed -fo ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_26.fa # bedtools: to FASTA format

	        ${NGStoolkitPath}/fa2kmerAbundanceTable.py -i ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_26.fa -k 2 -o ${mainPath}/simulation/XRseq/$1/control/$1_sorted_26_dinucleotideTable.txt # dinucleotide content of sequences at length 26

        fi

        if ${Key_sep_plus_minus}; then
        
	        grep "chr" ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted.bed | grep -v -e "chrY" -e "chrM" > ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_chr.bed

            mappedReads="$(grep -c '^' ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_chr.bed)"

            echo "Total Mapped Reads: ${mappedReads}" >> ${mainPath}/simulation/XRseq/$1/control/$1_control.txt # count

	        awk '{if($6=="+"){print}}' ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_chr.bed > ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_plus.bed # separating plus strand

	        awk '{if($6=="-"){print}}' ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_chr.bed > ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_minus.bed # separating minus strand

            echo bed files are ready!

        fi

        if ${Key_bedgraph_BigWig}; then
        
	        bedtools genomecov -i ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_plus.bed -g ${genomePath}/genome.fa.fai -bg -scale $(echo ${mappedReads} | awk '{print 1000000/$1}') > ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_plus.bdg # bedtools: to generate bedgraph files

	        bedtools genomecov -i ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_minus.bed -g ${genomePath}/genome.fa.fai -bg -scale $(echo ${mappedReads} | awk '{print -1000000/$1}') > ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_minus.bdg # bedtools: to generate bedgraph files

	        bedGraphToBigWig ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_plus.bdg ${genomePath}/genome.fa.fai ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_plus.bw # bedgraph to BigWig

	        bedGraphToBigWig ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_minus.bdg ${genomePath}/genome.fa.fai ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_minus.bw # bedgraph to BigWig

            echo BigWig files are ready!

        fi

        if ${Key_TS_NTS}; then

            bedtools intersect -sorted -a ${genomePath}/ensembl_genes.bed -b ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_chr.bed -wa -c -S -F 0.5 > ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_TScount.txt

            bedtools intersect -sorted -a ${genomePath}/ensembl_genes.bed -b ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_chr.bed -wa -c -s -F 0.5 > ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_NTScount.txt

            paste ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_TScount.txt ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_NTScount.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$14}' > ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_TSoverNTScount.txt

            echo TS/NTS counts are calculated!
        
        fi
    
    fi


# Analysis

    if ${Key_downstream_analysis}; then

        for ((i=0;i<${#dataset[@]};i++)); do

            echo data: ${dataset[i]} long name: ${data_name[i]} combine options: ${combine_options[i]}

            if ${Key_alignment}; then


                bedtools intersect -a ${mainPath}/data/${data_name[i]} -b ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_plus.bed -wa -c -F 0.5 >  ${mainPath}/simulation/XRseq/$1/intersect_combine/$1_sorted_plus_${dataset[i]}.txt

                bedtools intersect -a ${mainPath}/data/${data_name[i]} -b ${mainPath}/simulation/XRseq/$1/pre_analysis/$1_sorted_minus.bed -wa -c -F 0.5 >  ${mainPath}/simulation/XRseq/$1/intersect_combine/$1_sorted_minus_${dataset[i]}.txt

                echo ${dataset[i]} alignment done!
        
            fi

            if [[ ${dataset[i]} == *"windows"* ]]; then
                
                if ${Key_combineWindows}; then

		            ${mainPath}/combinewindows.py -i ${mainPath}/simulation/XRseq/$1/intersect_combine/$1_sorted_plus_${dataset[i]}.txt ${combine_options[i]} -o ${mainPath}/simulation/XRseq/$1/intersect_combine/$1_sorted_plus_${dataset[i]}_combined.txt

		            ${mainPath}/combinewindows.py -i ${mainPath}/simulation/XRseq/$1/intersect_combine/$1_sorted_minus_${dataset[i]}.txt ${combine_options[i]} -o ${mainPath}/simulation/XRseq/$1/intersect_combine/$1_sorted_minus_${dataset[i]}_combined.txt

                    echo ${dataset[i]} combined!

                fi
            
                if ${Key_moreInfo}; then

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/simulation/XRseq/$1/intersect_combine/$1_sorted_plus_${dataset[i]}_combined.txt -o ${mainPath}/simulation/XRseq/$1/final/$1_sorted_plus_${dataset[i]}_full.txt -c "${moreinfo}_sim" + "${mappedReads}"

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/simulation/XRseq/$1/intersect_combine/$1_sorted_minus_${dataset[i]}_combined.txt -o ${mainPath}/simulation/XRseq/$1/final/$1_sorted_minus_${dataset[i]}_full.txt -c "${moreinfo}_sim" - "${mappedReads}"

                    echo ${dataset[i]} info added!

                fi

            else 

                if ${Key_moreInfo}; then

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/simulation/XRseq/$1/intersect_combine/$1_sorted_plus_${dataset[i]}.txt -o ${mainPath}/simulation/XRseq/$1/final/$1_sorted_plus_${dataset[i]}_full.txt -c "${moreinfo}_sim" + "${mappedReads}"

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/simulation/XRseq/$1/intersect_combine/$1_sorted_minus_${dataset[i]}.txt -o ${mainPath}/simulation/XRseq/$1/final/$1_sorted_minus_${dataset[i]}_full.txt -c "${moreinfo}_sim" - "${mappedReads}"

                    echo ${dataset[i]} info added!

                fi

            fi

            if ${Key_rpkm}; then
            
                python ${mainPath}/RPKM.py -i ${mainPath}/simulation/XRseq/$1/final/$1_sorted_plus_${dataset[i]}_full.txt -chse 2 3 -c 7 -mr 0 -o ${mainPath}/simulation/XRseq/$1/final/$1_sorted_plus_${dataset[i]}_rpkm.txt

                python ${mainPath}/RPKM.py -i ${mainPath}/simulation/XRseq/$1/final/$1_sorted_minus_${dataset[i]}_full.txt -chse 2 3 -c 7 -mr 0 -o ${mainPath}/simulation/XRseq/$1/final/$1_sorted_minus_${dataset[i]}_rpkm.txt

                echo ${dataset[i]} rpkm calculated!

            fi

            if ${Key_file}; then

                cat ${mainPath}/simulation/XRseq/$1/final/$1_sorted_plus_${dataset[i]}_rpkm.txt >> ${mainPath}/${now}final_report_${dataset[i]}_sim.txt

                cat ${mainPath}/simulation/XRseq/$1/final/$1_sorted_minus_${dataset[i]}_rpkm.txt >> ${mainPath}/${now}final_report_${dataset[i]}_sim.txt 

                echo ${dataset[i]} final report created!

            fi

        done
        
    fi

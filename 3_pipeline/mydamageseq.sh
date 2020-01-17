#!/bin/bash

source ~/repairReplication/source_dir.sh
source ${mainPath}/source_key.sh
source ${mainPath}/source_dataset.sh
source ${mainPath}/functions_repairRep.sh

# Primary variables 

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

    mkdir -p ${mainPath}/Damageseq/$1/control # directory for control files
    mkdir -p ${mainPath}/Damageseq/$1/pre_analysis # directory for pre_analysis files
    mkdir -p ${mainPath}/Damageseq/$1/intersect_combine # directory for intersected and combined files
    mkdir -p ${mainPath}/Damageseq/$1/final # directory for final files

# Pre-analysis

    if ${Key_pre_analysis}; then

        check_raw_data $1 $rawdataPath  

        if [ $layout == "p" ]; then 

            TotalFastqReads_R1="$(grep -c '@' ${rawdataPath}/$1${zip})"

            TotalFastqReads_R2="$(grep -c '@' ${rawdataPath}/$1${zip})"

            TotalFastqReads="$(($TotalFastqReads_R1+$TotalFastqReads_R2))"

	        echo "Total Fastq Reads: ${TotalFastqReads}" > ${mainPath}/Damageseq/$1/control/$1_control.txt # Control: add total fastq reads     

            if ${Key_cutadapt}; then # cutadapt Paired     

                (cutadapt --discard-trimmed -g GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT -G GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT -o ${mainPath}/Damageseq/$1/pre_analysis/$1_R1_cutadapt.fastq.gz -p ${mainPath}/Damageseq/$1/pre_analysis/$1_R2_cutadapt.fastq.gz ${rawdataPath}/$1${zip} ${rawdataPath}/$1${zip}) >> ${mainPath}/Damageseq/$1/control/$1_control.txt

                TotalFilteredReads_R1="$(grep -c '@' ${mainPath}/Damageseq/$1/pre_analysis/$1_R1_cutadapt.fastq.gz)"

                TotalFilteredReads_R2="$(grep -c '@' ${mainPath}/Damageseq/$1/pre_analysis/$1_R2_cutadapt.fastq.gz)"

                TotalFilteredReads="$(($TotalFilteredReads_R1+$TotalFilteredReads_R2))"

	            echo "Total Reads after cutadapt filtering: ${TotalFilteredReads}" >> ${mainPath}/Damageseq/$1/control/$1_control.txt 

                echo cutadapt is done!

            fi

            if ${Key_bowtie2}; then # bowtie2 Paired
	
	            (bowtie2 -p 4 -X 1000 -x ${genomePath}/Bowtie2/genome -1 ${mainPath}/Damageseq/$1/pre_analysis/$1_R1_cutadapt.fastq.gz -2 ${mainPath}/Damageseq/$1/pre_analysis/$1_R2_cutadapt.fastq.gz -S ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.sam) 2>> ${mainPath}/Damageseq/$1/control/$1_control.txt

	            samtools view -q 20 -b -o ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.bam ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.sam # samtools: sam to bam

	            bedtools bamtobed -i ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.bam > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.bed # bedtools: bam to bed

                rm ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.sam
                rm ${mainPath}/Damageseq/$1/pre_analysis/$1_R1_cutadapt.fastq.gz
                rm ${mainPath}/Damageseq/$1/pre_analysis/$1_R2_cutadapt.fastq.gz

                echo bowtie2 mapping is done!

            fi

        elif [ $layout == "s" ]; then

            TotalFastqReads="$(grep -c '@' ${rawdataPath}/$1${zip})"
	
	        echo "Total Fastq Reads: ${TotalFastqReads}" > ${mainPath}/Damageseq/$1/control/$1_control.txt # Control: add total fastq reads 

            if ${Key_cutadapt}; then # cutadapt Single

                (cutadapt --discard-trimmed -g GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT -o ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.fastq.gz ${rawdataPath}/$1${zip}) >>  ${mainPath}/Damageseq/$1/control/$1_control.txt 

	            TotalFilteredReads="$(grep -c '@' ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.fastq.gz)"

	            echo "Total Reads after cutadapt filtering: ${TotalFilteredReads}" >> ${mainPath}/Damageseq/$1/control/$1_control.txt

                echo cutadapt is done! 

            fi

            if ${Key_bowtie2}; then
        
                (bowtie2 -p 4 -x ${genomePath}/Bowtie2/genome -U ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.fastq.gz -S ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.sam) 2>> ${mainPath}/Damageseq/$1/control/$1_control.txt

	            samtools view -q 20 -b -o ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.bam ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.sam # samtools: sam to bam

	            bedtools bamtobed -i ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.bam > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.bed # bedtools: bam to bed

                rm ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.sam
                rm ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.fastq.gz

                echo bowtie2 mapping is done!

            fi

        fi

        if ${Key_sort_count}; then

	        sort -u -k1,1 -k2,2n -k3,3n ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt.bed > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted.bed # sort


        fi

        if ${Key_sep_plus_minus}; then

	        grep "chr" ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted.bed | grep -v -e "chrY" -e "chrM" > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed

	        awk '{if($6=="+"){print}}' ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus.bed # separating plus strand

	        awk '{if($6=="-"){print}}' ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus.bed # separating minus strand

            bedtools flank -i ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus.bed -g ${genomePath}/genome.fa.fai -l 6 -r 0 > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_flanked_plus.bed # bedtools flank (each strand is 50bp)

	        bedtools flank -i ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus.bed -g ${genomePath}/genome.fa.fai -l 0 -r 6 > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_flanked_minus.bed # bedtools flank (each strand is 50bp)

	        bedtools slop -i ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_flanked_plus.bed -g ${genomePath}/genome.fa.fai -l 0 -r 4 > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_slopped_plus.bed # bedtools slop

	        bedtools slop -i ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_flanked_minus.bed -g ${genomePath}/genome.fa.fai -l 4 -r 0 > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_slopped_minus.bed # bedtools slop

	        awk '{print $3-$2}' ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_slopped_plus.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}' > ${mainPath}/Damageseq/$1/control/$1_cutadapt_plus_length_distribution.txt # length distribution

	        awk '{print $3-$2}' ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_slopped_minus.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}' > ${mainPath}/Damageseq/$1/control/$1_cutadapt_minus_length_distribution.txt # length distribution

	        awk '{ if ($3-$2 == 10) { print } }' ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_slopped_plus.bed > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus_10.bed # get only 10 nucleotide long

	        awk '{ if ($3-$2 == 10) { print } }' ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_slopped_minus.bed > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus_10.bed # get only 10 nucleotide long

	        bedtools getfasta -fi ${genomePath}/genome.fa -bed ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus_10.bed -fo ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus_10.fa -s # bedtools: to FASTA format

	        bedtools getfasta -fi ${genomePath}/genome.fa -bed ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus_10.bed -fo ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus_10.fa -s # bedtools: to FASTA format

	        ${NGStoolkitPath}/fa2kmerAbundanceTable.py -i ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus_10.fa -k 2 -o ${mainPath}/Damageseq/$1/control/$1_cutadapt_sorted_plus_10_dinucleotideTable.txt # dinucleotide content

	        ${NGStoolkitPath}/fa2kmerAbundanceTable.py -i ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus_10.fa -k 2 -o ${mainPath}/Damageseq/$1/control/$1_cutadapt_sorted_minus_10_dinucleotideTable.txt # dinucleotide content

	        awk '{print $1"\t"$6}' ${mainPath}/Damageseq/$1/control/$1_cutadapt_sorted_plus_10_dinucleotideTable.txt > ${mainPath}/Damageseq/$1/control/$1_cutadapt_sorted_plus_10_dinucleotideTable_pos5-6.txt # damage positions

	        awk '{print $1"\t"$6}' ${mainPath}/Damageseq/$1/control/$1_cutadapt_sorted_minus_10_dinucleotideTable.txt > ${mainPath}/Damageseq/$1/control/$1_cutadapt_sorted_minus_10_dinucleotideTable_pos5-6.txt # damage positions

	        ${NGStoolkitPath}/fa2bedByChoosingReadMotifs.py -i ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus_10.fa -o ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus_dipyrimidines.bed -r ".{4}(c|t|C|T){2}.{4}" # taking only dipyrimidines

	        ${NGStoolkitPath}/fa2bedByChoosingReadMotifs.py -i ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus_10.fa -o ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus_dipyrimidines.bed -r ".{4}(c|t|C|T){2}.{4}" # taking only dipyrimidines
        
            minus_line="$(grep -c '^' ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus_dipyrimidines.bed)"
            plus_line="$(grep -c '^' ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus_dipyrimidines.bed)"
            mappedReads=`echo "$minus_line + $plus_line" | bc`

            echo "Total Mapped Reads: ${mappedReads}" >> ${mainPath}/Damageseq/$1/control/$1_control.txt # count
    
            echo bed files are ready!

        fi

        if ${Key_bedgraph_BigWig}; then     

	        bedtools genomecov -i ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus.bed -g ${genomePath}/genome.fa.fai -bg -scale $(echo ${mappedReads} | awk '{print 1000000/$1}') > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus.bdg # bedtools: to generate bedgraph files

	        bedtools genomecov -i ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus.bed -g ${genomePath}/genome.fa.fai -bg -scale $(echo ${mappedReads} | awk '{print -1000000/$1}') > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus.bdg # bedtools: to generate bedgraph files

	        bedGraphToBigWig ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus.bdg ${genomePath}/genome.fa.fai ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus.bw # bedgraph to BigWig

	        bedGraphToBigWig ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus.bdg ${genomePath}/genome.fa.fai ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus.bw # bedgraph to BigWig

            echo BigWig files are ready!

        fi

        if ${Key_TS_NTS}; then

            cat ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus_dipyrimidines.bed ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus_dipyrimidines.bed | sort -k1,1 -k2,2n -k3,3n > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_dipyrimidines.bed

            bedtools intersect -sorted -a ${genomePath}/ensembl_genes.bed -b ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_dipyrimidines.bed -wa -c -S -F 0.5 > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_TScount.txt

            bedtools intersect -sorted -a ${genomePath}/ensembl_genes.bed -b ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_dipyrimidines.bed -wa -c -s -F 0.5 > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_NTScount.txt

            paste ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_TScount.txt ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_NTScount.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$14}' > ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_TSoverNTScount.txt         

            echo TS/NTS counts are calculated!   
        
        fi
    fi

# Analysis

    if ${Key_downstream_analysis}; then

        for ((i=0;i<${#dataset[@]};i++)); do

            echo data: ${dataset[i]} long name: ${data_name[i]} combine options: ${combine_options[i]}

            if ${Key_alignment}; then


                bedtools intersect -a ${mainPath}/data/${data_name[i]} -b ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_plus_dipyrimidines.bed -wa -c -F 0.5 >  ${mainPath}/Damageseq/$1/intersect_combine/$1_cutadapt_sorted_plus_${dataset[i]}.txt

                bedtools intersect -a ${mainPath}/data/${data_name[i]} -b ${mainPath}/Damageseq/$1/pre_analysis/$1_cutadapt_sorted_minus_dipyrimidines.bed -wa -c -F 0.5 >  ${mainPath}/Damageseq/$1/intersect_combine/$1_cutadapt_sorted_minus_${dataset[i]}.txt

                echo ${dataset[i]} alignment done!
        
            fi

            if [[ ${dataset[i]} == *"windows"* ]]; then
                
                if ${Key_combineWindows}; then

		            ${mainPath}/combinewindows.py -i ${mainPath}/Damageseq/$1/intersect_combine/$1_cutadapt_sorted_plus_${dataset[i]}.txt ${combine_options[i]} -o ${mainPath}/Damageseq/$1/intersect_combine/$1_cutadapt_sorted_plus_${dataset[i]}_combined.txt

		            ${mainPath}/combinewindows.py -i ${mainPath}/Damageseq/$1/intersect_combine/$1_cutadapt_sorted_minus_${dataset[i]}.txt ${combine_options[i]} -o ${mainPath}/Damageseq/$1/intersect_combine/$1_cutadapt_sorted_minus_${dataset[i]}_combined.txt

                    echo ${dataset[i]} combined!

                fi
            
                if ${Key_moreInfo}; then

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/Damageseq/$1/intersect_combine/$1_cutadapt_sorted_plus_${dataset[i]}_combined.txt -o ${mainPath}/Damageseq/$1/final/$1_cutadapt_sorted_plus_${dataset[i]}_full.txt -c "${moreinfo}" + "${mappedReads}"

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/Damageseq/$1/intersect_combine/$1_cutadapt_sorted_minus_${dataset[i]}_combined.txt -o ${mainPath}/Damageseq/$1/final/$1_cutadapt_sorted_minus_${dataset[i]}_full.txt -c "${moreinfo}" - "${mappedReads}"

                    echo ${dataset[i]} info added!

                fi

            else 

                if ${Key_moreInfo}; then

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/Damageseq/$1/intersect_combine/$1_cutadapt_sorted_plus_${dataset[i]}.txt -o ${mainPath}/Damageseq/$1/final/$1_cutadapt_sorted_plus_${dataset[i]}_full.txt -c "${moreinfo}" + "${mappedReads}"

                    ${NGStoolkitPath}/addColumns.py -i ${mainPath}/Damageseq/$1/intersect_combine/$1_cutadapt_sorted_minus_${dataset[i]}.txt -o ${mainPath}/Damageseq/$1/final/$1_cutadapt_sorted_minus_${dataset[i]}_full.txt -c "${moreinfo}" - "${mappedReads}"

                    echo ${dataset[i]} info added!

                fi

            fi

            if ${Key_rpkm}; then
            
                python ${mainPath}/RPKM.py -i ${mainPath}/Damageseq/$1/final/$1_cutadapt_sorted_plus_${dataset[i]}_full.txt -chse 2 3 -c 7 -mr 0 -o ${mainPath}/Damageseq/$1/final/$1_cutadapt_sorted_plus_${dataset[i]}_rpkm.txt

                python ${mainPath}/RPKM.py -i ${mainPath}/Damageseq/$1/final/$1_cutadapt_sorted_minus_${dataset[i]}_full.txt -chse 2 3 -c 7 -mr 0 -o ${mainPath}/Damageseq/$1/final/$1_cutadapt_sorted_minus_${dataset[i]}_rpkm.txt

                echo ${dataset[i]} rpkm calculated!

            fi

            if ${Key_file}; then

                cat ${mainPath}/Damageseq/$1/final/$1_cutadapt_sorted_plus_${dataset[i]}_rpkm.txt >> ${mainPath}/${now}final_report_${dataset[i]}.txt

                cat ${mainPath}/Damageseq/$1/final/$1_cutadapt_sorted_minus_${dataset[i]}_rpkm.txt >> ${mainPath}/${now}final_report_${dataset[i]}.txt 

                echo ${dataset[i]} final report created!

            fi

        done
        
    fi

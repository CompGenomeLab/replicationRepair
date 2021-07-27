#!/bin/bash

# This function checks the path, layout and compression (gz zipped or not) of the file.
check_raw_data() {

    local sample=$1
    local rawPath=$2

    if [[ $(find "${rawPath}" -name "${sample}_R1.fastq.gz" -o -name "${sample}_1.fastq.gz" ) ]]; then
        rawdataPath=$(find "${rawPath}" -name "${sample}_R1.fastq.gz" -o -name "${sample}_1.fastq.gz" | sed -ne "s/ *\/${sample}.*//p" | head -1)
        zip=$(find "${rawPath}" -name "${sample}_R1.fastq.gz" -o -name "${sample}_1.fastq.gz" | sed -ne "s/ *.*${sample}//p" | head -1)
        zip2=$(echo ${zip} | tr 1 2)
        layout="p"
        echo raw data is compressed! layout is paired!
    elif [[ $(find "${rawPath}" -name "${sample}.fastq.gz") ]]; then
        rawdataPath=$(find "${rawPath}" -name "${sample}.fastq.gz" | sed -ne "s/ *\/${sample}.*//p" | head -1) 
        zip=$(find "${rawPath}" -name "${sample}.fastq.gz" | sed -ne "s/ *.*${sample}//p" | head -1)
        layout="s"
        echo raw data is compressed! layout is single!
    elif [[ $(find "${rawPath}" -name "${sample}_R1.fastq" -o -name "${sample}_1.fastq") ]]; then 
        rawdataPath=$(find "${rawPath}" -name "${sample}_R1.fastq" -o -name "${sample}_1.fastq" | sed -ne "s/ *\/${sample}.*//p" | head -1)
        zip=$(find "${rawPath}" -name "${sample}_R1.fastq" -o -name "${sample}_1.fastq" | sed -ne "s/ *.*${sample}//p" | head -1)
        zip2=$(echo ${zip} | tr 1 2)
        layout="p"
        echo raw data is not compressed! layout is paired!   
    elif [[ $(find "${rawPath}" -name "${sample}.fastq") ]]; then 
        rawdataPath=$(find "${rawPath}" -name "${sample}.fastq" | sed -ne "s/ *\/${sample}.*//p" | head -1)
        zip=$(find "${rawPath}" -name "${sample}.fastq" | sed -ne "s/ *.*${sample}//p" | head -1)
        layout="s"
        echo raw data is not compressed! layout is single! 
    else      
        echo raw data is missing!
    fi
}


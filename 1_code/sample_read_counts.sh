#!/bin/bash

source ~/repairReplication/source_dir.sh

samples=(NT4H1-CT5H2-HelaD1-5R2h2_combined XT4H1-CT5H2-HelaD3-5R2h2_combined R19029847-HD2TCD4-HelaD1-5R2h1_combined R19029847-HD2TCD4-HelaD3-5R2h1_combined HDA64A1_ATCACG HDA64B19_GTGAAA HDE64A4_TGACCA HDE64B20_GTGGCC HDL64A5_ACAGTG HDL64B22_CGTACG HDACA6_GCCAAT HDACB23_GAGTGG HDECA10_TAGCTT HDECB25_ACTGAT HDLCA12_CTTGTA HDLCB27_ATTCCT HXA64A1_ATCACG HXA64B7_CAGATC HXE64A2_CGATGT HXE64B8_ACTTGA HXL64A3_TTAGGC HXL64B9_GATCAG HXACA4_TGACCA HXACB10_TAGCTT HXECA5_ACAGTG HXECB11_GGCTAC HXLCA6_GCCAAT HXLCB12_CTTGTA R19026421-2019XR1-Hela15X2_combined_R1 R19026421-2019XR1-Hela35X3_combined_R1 R19033030-2019XR3-Hela15X7_combined_R1 R19033030-2019XR3-Hela35X8_combined_R1 R19040024-HelainputD-Hl15R2hID_combined R19040024-HelainputD-Hl0h20JID_combined R19040024-HelainputD-Hl35R2hID_combined)

echo "sample_name,read" > ${mainPath}/read_counts.csv

for ((i=0;i<${#samples[@]};i++)); do

    getname="$(echo ${samples[i]} | sed 's/-/_/g')"
    moreinfo="$(grep ${getname} ${mainPath}/project_repair_replication_all_samples.csv | sed 's/,/\t/g' | awk '{print $1}')"

    if [ -f ${mainPath}/Input/${samples[i]}/pre_analysis/${samples[i]}_chr.bed ]; then
	    mappedReads="$(grep -c '^' ${mainPath}/Input/${samples[i]}/pre_analysis/${samples[i]}_chr.bed)"

    elif [ -f ${mainPath}/Damageseq/${samples[i]}/pre_analysis/${samples[i]}_cutadapt_sorted_minus_dipyrimidines.bed ]; then
        minus_line="$(grep -c '^' ${mainPath}/Damageseq/${samples[i]}/pre_analysis/${samples[i]}_cutadapt_sorted_minus_dipyrimidines.bed)"
        plus_line="$(grep -c '^' ${mainPath}/Damageseq/${samples[i]}/pre_analysis/${samples[i]}_cutadapt_sorted_plus_dipyrimidines.bed)"
        mappedReads=`echo "$minus_line + $plus_line" | bc`

    elif [ -f ${mainPath}/XRseq/${samples[i]}/pre_analysis/${samples[i]}_cutadapt_sorted_chr.bed ]; then
	    mappedReads="$(grep -c '^' ${mainPath}/XRseq/${samples[i]}/pre_analysis/${samples[i]}_cutadapt_sorted_chr.bed)" 

    fi

    echo ${moreinfo},${mappedReads} >> ${mainPath}/read_counts.csv

done

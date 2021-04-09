awk '{print $0"\t""Termination_Zones""\t"".""\t""."}' GM_hmm_HMMsegments_TZ.bed | sort -k1,1 -k2,2n -k3,3n > GM_hmm_HMMsegments_TZ_org_sorted.bed

awk '{print $0"\t"".""\t""."}' GSE53984_GSM923443_Gm06990_Rep1_segments.bed | sort -k1,1 -k2,2n -k3,3n > GSE53984_GSM923443_Gm06990_Rep1_segments_org_sorted.bed

bedtools intersect -sorted -a GM_hmm_HMMsegments_TZ_org_sorted.bed -b GSE53984_GSM923443_Gm06990_Rep1_segments_org_sorted.bed -wa -wb -f 0.5 | awk '{print $1"\t"$2"\t"$3"\t"$4"_"$10"\t"$5"\t"$6}' > TZ.gm06990_repdomains.bed

# for damage-seq
cat ${preAnalysis}/$1_cutadapt_sorted_minus_dipyrimidines.bed ${preAnalysis}/$1_cutadapt_sorted_plus_dipyrimidines.bed | sort -k1,1 -k2,2n -k3,3n > ${preAnalysis}/$1_cutadapt_sorted_dipyrimidines.bed

bedtools intersect -sorted -a ${mainPath}/data/hg19_ucsc_genes_knownCanonical_tss_windows_201_100.bed -b ${preAnalysis}/$1_cutadapt_sorted_dipyrimidines.bed -wa -c -S -F 0.5 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""TS""\t"$7}' > ${preAnalysis}/$1_cutadapt_sorted_TScount.txt

bedtools intersect -sorted -a ${mainPath}/data/hg19_ucsc_genes_knownCanonical_tss_windows_201_100.bed -b ${preAnalysis}/$1_cutadapt_sorted_dipyrimidines.bed -wa -c -s -F 0.5 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""NTS""\t"$7}' > ${preAnalysis}/$1_cutadapt_sorted_NTScount.txt

cat ${preAnalysis}/$1_cutadapt_sorted_TScount.txt ${preAnalysis}/$1_cutadapt_sorted_NTScount.txt > ${preAnalysis}/$1_cutadapt_sorted_TS_NTScount.txt

${mainPath}/combinewindows.py -i ${preAnalysis}/$1_cutadapt_sorted_TS_NTScount.txt -strand T -o ${preAnalysis}/$1_cutadapt_sorted_TS_NTScount_combined.txt

${NGStoolkitPath}/addColumns.py -i ${preAnalysis}/$1_cutadapt_sorted_TS_NTScount_combined.txt -o ${preAnalysis}/$1_cutadapt_sorted_TS_NTScount_combined_full.txt -c "${moreinfo}" "." "${mappedReads}"

python ${mainPath}/RPKM.py -i ${preAnalysis}/$1_cutadapt_sorted_TS_NTScount_combined_full.txt -chse 2 3 -c 7 -mr 0 -o ${preAnalysis}/$1_cutadapt_sorted_TS_NTScount_combined_rpkm.txt

cat ${preAnalysis}/$1_cutadapt_sorted_TS_NTScount_combined_rpkm.txt >> ${mainPath}/final_report_tss.txt


# for xr-seq
bedtools intersect -sorted -a ${mainPath}/data/hg19_ucsc_genes_knownCanonical_tss_windows_201_100.bed -b ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed -wa -c -S -F 0.5 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""TS""\t"$7}' > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TScount.txt

bedtools intersect -sorted -a ${mainPath}/data/hg19_ucsc_genes_knownCanonical_tss_windows_201_100.bed -b ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_chr.bed -wa -c -s -F 0.5 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""NTS""\t"$7}' > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_NTScount.txt

cat ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TScount.txt ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_NTScount.txt > ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TS_NTScount.txt

${mainPath}/combinewindows.py -i ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TS_NTScount.txt -strand T -o ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TS_NTScount_combined.txt

${NGStoolkitPath}/addColumns.py -i ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TS_NTScount_combined.txt -o ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TS_NTScount_combined_full.txt -c "${moreinfo}" "." "${mappedReads}"

python ${mainPath}/RPKM.py -i ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TS_NTScount_combined_full.txt -chse 2 3 -c 7 -mr 0 -o ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TS_NTScount_combined_rpkm.txt

cat ${mainPath}/XRseq/$1/pre_analysis/$1_cutadapt_sorted_TS_NTScount_combined_rpkm.txt >> ${mainPath}/final_report_tss.txt
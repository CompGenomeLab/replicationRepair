# set directories
outputPath="/home/azgarian/Desktop/ssm"
filterPath="/home/azgarian/Desktop/masaustu"
genome="/home/azgarian/Documents/myprojects/\
replicationRepair/0_data/gitignore/Genome/\
hg19_genome.fa"
mutationFile="/home/azgarian/Documents/myprojects/replicationRepair/0_data/\
gitignore/Original/simple_somatic_mutation.open.MELA-AU.tsv.gz"
newname="melanoma"
toolPath="/home/azgarian/Documents/myprojects/replicationRepair/1_code/python"
dataPath="/home/azgarian/Documents/myprojects/replicationRepair/0_data/gitignore"
finalPath="/home/azgarian/Documents/myprojects/replicationRepair/3_output/gitignore/1_TextforPlotting"


# unzipping .gz file
#echo `date +%R:%S`: unzipping..
#zcat ${mutationFile} | \
#    sed 's/single base sub/single_base_sub/g' > \
#    ${outputPath}/${newname}.tsv
#echo `date +%R:%S`: done!


# extracting single base substitution mutations
#echo `date +%R:%S`: get only subs..
#awk -v t="\t" '{print $2t"chr"$9t$10t$11t$14t$16t$17}' \
#    ${outputPath}/${newname}.tsv | \
#    sort -k1,1 -k2,2 -k3,3n -k4,4n | \
#    uniq | \
#    grep single_base_sub > \
#    ${outputPath}/${newname}_only_subs.tsv
#echo `date +%R:%S`: done!


# organize format
#echo `date +%R:%S`: organize..
#grep -v -e "chrY" -e "chrMT" ${outputPath}/${newname}_only_subs.tsv | \
#    awk '{print $2"\t"($3-2)"\t"($4+1)"\t"$6"_"$7"\t""0""\t""+"}' > \
#    ${outputPath}/${newname}_only_subs_org.tsv
#echo `date +%R:%S`: done!


# create fasta file
#echo `date +%R:%S`: creating fasta file..
#bedtools getfasta \
#    -fi ${genome} \
#    -bed ${outputPath}/${newname}_only_subs_org.tsv \
#    -fo ${outputPath}/${newname}.fa \
#    -name 
#echo `date +%R:%S`: done!


# create desired .tsv format, sorted  
#echo `date +%R:%S`: creating final bed format..
#awk -F'[>_:-]' \
#    'NR%2==1{a=$5"\t"$6"\t"$7"\t"$2"_"$3} NR%2==0{print a"_"$0}' \
#    ${outputPath}/${newname}.fa | \
#    sort -k1,1 -k2,2n -k3,3n > \
#    ${outputPath}/${newname}.bed
#echo `date +%R:%S`: done!


# filter C->T and G->A mutations
#echo `date +%R:%S`: filtering target mutations..
#grep -e "C_T_[TC]" -e "G_A_..[GA]" \
#    ${outputPath}/${newname}.bed | \
#    awk -v t="\t" '{\
#    if ($4 ~ "C_T"){print $0t"0"t"+"}; \
#    if ($4 ~ "G_A"){print $0t"0"t"-"}}' > \
#    ${outputPath}/${newname}_target_mut.tsv
#echo `date +%R:%S`: done!

# separate plus and minus
#echo `date +%R:%S`: separate plus and minus..
#awk '{if($6=="+"){print}}' ${outputPath}/${newname}_target_mut.tsv > \
#    ${outputPath}/${newname}_target_mut_plus.tsv
#awk '{if($6=="-"){print}}' ${outputPath}/${newname}_target_mut.tsv > \
#    ${outputPath}/${newname}_target_mut_minus.tsv
#echo `date +%R:%S`: done!

## Example:
# convert to bed format
#echo `date +%R:%S`: converting to bed..
#awk -v t="\t" '{\
#    if ($4 ~ "TCG"){print $1t$2t$3t"."t"0"t"+"}; \
#    if ($4 ~ "CGA"){print $1t$2t$3t"."t"0"t"-"}}' \
#    ${outputPath}/${newname}_target_mut.tsv > \
#    ${outputPath}/${newname}_target_mut_TCG.bed
#echo `date +%R:%S`: done!

#echo `date +%R:%S`: converting to bed..
#awk -v t="\t" '{\
#    if ($4 ~ "CCG"){print $1t$2t$3t"."t"0"t"+"}; \
#    if ($4 ~ "CGG"){print $1t$2t$3t"."t"0"t"-"}}' \
#    ${outputPath}/${newname}_target_mut.tsv > \
#    ${outputPath}/${newname}_target_mut_CCG.bed
#echo `date +%R:%S`: done!

#echo `date +%R:%S`: converting to bed..
#awk -v t="\t" '{\
#    if ($4 ~ "CCA"){print $1t$2t$3t"."t"0"t"+"}; \
#    if ($4 ~ "TGG"){print $1t$2t$3t"."t"0"t"-"}}' \
#    ${outputPath}/${newname}_target_mut.tsv > \
#    ${outputPath}/${newname}_target_mut_CCA.bed
#echo `date +%R:%S`: done!

#echo `date +%R:%S`: converting to bed..
#awk -v t="\t" '{\
#    if ($4 ~ "CCT"){print $1t$2t$3t"."t"0"t"+"}; \
#    if ($4 ~ "AGG"){print $1t$2t$3t"."t"0"t"-"}}' \
#    ${outputPath}/${newname}_target_mut.tsv > \
#    ${outputPath}/${newname}_target_mut_CCT.bed
#echo `date +%R:%S`: done!

#echo `date +%R:%S`: converting to bed..
#awk -v t="\t" '{\
#    if ($4 ~ "CCC"){print $1t$2t$3t"."t"0"t"+"}; \
#    if ($4 ~ "GGG"){print $1t$2t$3t"."t"0"t"-"}}' \
#    ${outputPath}/${newname}_target_mut.tsv > \
#    ${outputPath}/${newname}_target_mut_CCC.bed
#echo `date +%R:%S`: done!

#echo `date +%R:%S`: converting to bed..
#awk -v t="\t" '{\
#    if ($4 ~ "TCT"){print $1t$2t$3t"."t"0"t"+"}; \
#    if ($4 ~ "AGA"){print $1t$2t$3t"."t"0"t"-"}}' \
#    ${outputPath}/${newname}_target_mut.tsv > \
#    ${outputPath}/${newname}_target_mut_TCT.bed
#echo `date +%R:%S`: done!

#echo `date +%R:%S`: converting to bed..
#awk -v t="\t" '{\
#    if ($4 ~ "TCA"){print $1t$2t$3t"."t"0"t"+"}; \
#    if ($4 ~ "TGA"){print $1t$2t$3t"."t"0"t"-"}}' \
#    ${outputPath}/${newname}_target_mut.tsv > \
#    ${outputPath}/${newname}_target_mut_TCA.bed
#echo `date +%R:%S`: done!

#echo `date +%R:%S`: converting to bed..
#awk -v t="\t" '{\
#    if ($4 ~ "TCC"){print $1t$2t$3t"."t"0"t"+"}; \
#    if ($4 ~ "GGA"){print $1t$2t$3t"."t"0"t"-"}}' \
#    ${outputPath}/${newname}_target_mut.tsv > \
#    ${outputPath}/${newname}_target_mut_TCC.bed
#echo `date +%R:%S`: done!


################################################################################

now=$(date +"[%Y.%m.%d.%H:%M]") # date

bedtools intersect -a ${outputPath}/${newname}_target_mut_plus.tsv -b ${filterPath}/hg19_ucsc_genes_org.bed -v -f 0.5 >  ${outputPath}/${newname}_target_mut_plus_filtered.tsv

bedtools intersect -a ${outputPath}/${newname}_target_mut_minus.tsv -b ${filterPath}/hg19_ucsc_genes_org.bed -v -f 0.5 >  ${outputPath}/${newname}_target_mut_minus_filtered.tsv

bedtools intersect -a ${dataPath}/InZones/initiation.zones.repdomains.hela_windows_201_100.bed -b ${outputPath}/${newname}_target_mut_plus_filtered.tsv -wa -c -F 0.5 > ${outputPath}/${newname}_target_mut_plus_intersect.bed

bedtools intersect -a ${dataPath}/InZones/initiation.zones.repdomains.hela_windows_201_100.bed -b ${outputPath}/${newname}_target_mut_minus_filtered.tsv -wa -c -F 0.5 > ${outputPath}/${newname}_target_mut_minus_intersect.bed

#bedtools intersect -a ${dataPath}/InZones/initiation.zones.repdomains.hela_windows_201_100.bed -b ${outputPath}/${newname}_target_mut_plus.tsv -wa -c -F 0.5 > ${outputPath}/${newname}_target_mut_plus_intersect.bed

#bedtools intersect -a ${dataPath}/InZones/initiation.zones.repdomains.hela_windows_201_100.bed -b ${outputPath}/${newname}_target_mut_minus.tsv -wa -c -F 0.5 > ${outputPath}/${newname}_target_mut_minus_intersect.bed

awk '{print $0"\t""+"}' ${outputPath}/${newname}_target_mut_plus_intersect.bed > ${outputPath}/${newname}_target_mut_combined.bed

awk '{print $0"\t""-"}' ${outputPath}/${newname}_target_mut_minus_intersect.bed >> ${outputPath}/${newname}_target_mut_combined.bed

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$8"\t"$7}' ${outputPath}/${newname}_target_mut_combined.bed > ${outputPath}/${newname}_target_mut_combined_org.bed

${toolPath}/combinewindows.py -i ${outputPath}/${newname}_target_mut_combined_org.bed -score T -strand T -o ${finalPath}/${now}inzones_repdomains_201_100_${newname}_mutations_combined_intergenic.bed
 



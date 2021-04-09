# set directories
source ../1_code/source_dir.sh

mutationFile="" # mutation file with full path 
data="" # regions data that mutations will be intersected to

mutPath="${mainPath}/0_data/Melanoma"

# unzipping .gz file
echo `date +%R:%S`: unzipping..
zcat ${mutationFile} | \
    sed 's/single base sub/single_base_sub/g' > \
    ${mutPath}/melanoma.tsv
echo `date +%R:%S`: done!

# extracting single base substitution mutations
echo `date +%R:%S`: get only subs..
awk -v t="\t" '{print $2t"chr"$9t$10t$11t$14t$16t$17}' \
    ${mutPath}/melanoma.tsv | \
    sort -k1,1 -k2,2 -k3,3n -k4,4n | \
    uniq | \
    grep single_base_sub > \
    ${mutPath}/melanoma_only_subs.tsv
echo `date +%R:%S`: done!

# organize format
echo `date +%R:%S`: organize..
grep -v -e "chrY" -e "chrMT" ${mutPath}/melanoma_only_subs.tsv | \
    awk '{print $2"\t"($3-2)"\t"($4+1)"\t"$6"_"$7"\t""0""\t""+"}' > \
    ${mutPath}/melanoma_only_subs_org.tsv
echo `date +%R:%S`: done!

# create fasta file
echo `date +%R:%S`: creating fasta file..
bedtools getfasta \
    -fi ${genomePath}/genome.fa \
    -bed ${mutPath}/melanoma_only_subs_org.tsv \
    -fo ${mutPath}/melanoma.fa \
    -name 
echo `date +%R:%S`: done!

# create desired .tsv format, sorted  
echo `date +%R:%S`: creating final bed format..
awk -F'[>_:-]' \
    'NR%2==1{a=$5"\t"$6"\t"$7"\t"$2"_"$3} NR%2==0{print a"_"$0}' \
    ${mutPath}/melanoma.fa | \
    sort -k1,1 -k2,2n -k3,3n > \
    ${mutPath}/melanoma.bed
echo `date +%R:%S`: done!

# filter C->T and G->A mutations
echo `date +%R:%S`: filtering target mutations..
grep -e "C_T_[TC]" -e "G_A_..[GA]" \
    ${mutPath}/melanoma.bed | \
    awk -v t="\t" '{\
    if ($4 ~ "C_T"){print $0t"0"t"+"}; \
    if ($4 ~ "G_A"){print $0t"0"t"-"}}' > \
    ${mutPath}/melanoma_target_mut.tsv
echo `date +%R:%S`: done!

# separate plus and minus
echo `date +%R:%S`: separate plus and minus..
awk '{if($6=="+"){print}}' ${mutPath}/melanoma_target_mut.tsv > \
    ${mutPath}/melanoma_target_mut_plus.tsv
awk '{if($6=="-"){print}}' ${mutPath}/melanoma_target_mut.tsv > \
    ${mutPath}/melanoma_target_mut_minus.tsv
echo `date +%R:%S`: done!

## If only a trinucleotide desired:
# convert to bed format
echo `date +%R:%S`: converting to bed..
awk -v t="\t" '{\
    if ($4 ~ "TCG"){print $1t$2t$3t"."t"0"t"+"}; \
    if ($4 ~ "CGA"){print $1t$2t$3t"."t"0"t"-"}}' \
    ${mutPath}/melanoma_target_mut.tsv > \
    ${mutPath}/melanoma_target_mut_TCG.bed
echo `date +%R:%S`: done!


################################################################################

now=$(date +"[%Y.%m.%d.%H:%M]") # date

# filtering genic mutations
bedtools intersect -a ${mutPath}/melanoma_target_mut_plus.tsv -b ${genomePath}/hg19_ucsc_genes_org.bed -v -f 0.5 >  ${mutPath}/melanoma_target_mut_plus_filtered.tsv

bedtools intersect -a ${mutPath}/melanoma_target_mut_minus.tsv -b ${genomePath}/hg19_ucsc_genes_org.bed -v -f 0.5 >  ${mutPath}/melanoma_target_mut_minus_filtered.tsv

# intersecting with the defined regions
bedtools intersect -a ${data} -b ${mutPath}/melanoma_target_mut_plus_filtered.tsv -wa -c -F 0.5 > ${mutPath}/melanoma_target_mut_plus_intersect.bed

bedtools intersect -a ${data} -b ${mutPath}/melanoma_target_mut_minus_filtered.tsv -wa -c -F 0.5 > ${mutPath}/melanoma_target_mut_minus_intersect.bed

# combine plus and minus strands
awk '{print $0"\t""+"}' ${mutPath}/melanoma_target_mut_plus_intersect.bed > ${mutPath}/melanoma_target_mut_combined.bed

awk '{print $0"\t""-"}' ${mutPath}/melanoma_target_mut_minus_intersect.bed >> ${mutPath}/melanoma_target_mut_combined.bed

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$8"\t"$7}' ${mutPath}/melanoma_target_mut_combined.bed > ${mutPath}/melanoma_target_mut_combined_org.bed

# if regions will be aggregated based on score, strand, and name
${codePath}/combinewindows.py -i ${mutPath}/melanoma_target_mut_combined_org.bed -score T -strand T -o ${mutPath}/${now}iz_hela_to_gm_imr_repdomains_with_scores_201_100_melanoma_mutations_combined_intergenic_agg.bed

# if no aggregation needed 
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$8"\t"$7}' ${mutPath}/melanoma_target_mut_combined.bed > ${mutPath}/${now}iz_hela_to_gm_201_100_melanoma_mutations_combined_intergenic.bed
 
#!/bin/bash
#
# CompecTA (c) 2018
#
# NAMD job submission script
#
# TODO:
#   - Set name of the job below changing "NAMD" value.
#   - Set the requested number of nodes (servers) with --nodes parameter.
#   - Set the requested number of tasks (cpu cores) with --ntasks parameter. (Total accross all nodes)
#   - Select the partition (queue) you want to run the job in:
#     - short : For jobs that have maximum run time of 120 mins. Has higher priority.
#     - mid   : For jobs that have maximum run time of 1 days. Lower priority than short.
#     - long  : For jobs that have maximum run time of 7 days. Lower priority than long.
#     - longer: For testing purposes, queue has 15 days limit but only 2 nodes.
#     - cuda  : For CUDA jobs. Solver that can utilize CUDA acceleration can use this queue. 7 days limit.
#   - Set the required time limit for the job with --time parameter.
#     - Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"
#   - Put this script and all the input file under the same directory.
#   - Set the required parameters, input/output file names below.
#   - If you do not want mail please remove the line that has --mail-type and --mail-user. If you do want to get notification emails, set your email address.
#   - Put this script and all the input file under the same directory.
#   - Submit this file using:
#      sbatch slurm_submit.sh
#
# -= Resources =-
#
#SBATCH --job-name=genome_windows
#SBATCH --account=users
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --qos=long
#SBATCH --partition=long
#SBATCH --time=7-00
#SBATCH --output=%j-slurm.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cemazgari@sabanciuniv.edu
#SBATCH --array=0-34


################################################################################
source /etc/profile.d/modules.sh
echo "source /etc/profile.d/modules.sh"
################################################################################

# Module File
#echo "Loading Foo..."
#module load foo-tksd6ij
#echo

echo ""
echo "======================================================================================"
env
echo "======================================================================================"
echo ""

echo "======================================================================================"
# Set stack size to unlimited
echo "Setting stack size to unlimited..."
ulimit -s unlimited
ulimit -l unlimited
ulimit -a
echo

echo "Running Example Job...!"
echo "==============================================================================="
# Command 1 for matrix
echo "Running Python script..."
# Put Python script command below

# Command 2 for matrix
echo "Running G++ compiler..."
# Put g++ compiler command below

# Command 3 for matrix
echo "Running compiled binary..."
# Put compiled binary command below

module load bowtie2-2.3.4.1-gcc-8.2.0-dqmpp4a
module load bedtools2-2.27.1-gcc-7.2.0-tn3naqh
module load bamtools-2.5.1-gcc-8.2.0-2quxjro
module load samtools-1.9-gcc-8.2.0-r22yn5w
module load bioconda-adebali
module load pandas-0.24.1



mainPath="/cta/users/cemazgari/repairReplication"
NGStoolkitPath="/cta/users/cemazgari/NGStoolkit/bin"

echo This is task $SLURM_ARRAY_TASK_ID


samples=(R19026421-2019XR1-Hela15X2_combined_R1 R19026421-2019XR1-Hela35X3_combined_R1 R19033030-2019XR3-Hela15X7_combined_R1 R19033030-2019XR3-Hela35X8_combined_R1 NT4H1-CT5H2-HelaD1-5R2h2_combined XT4H1-CT5H2-HelaD3-5R2h2_combined R19029847-HD2TCD4-HelaD1-5R2h1_combined R19029847-HD2TCD4-HelaD3-5R2h1_combined HXA64A1_ATCACG HXA64B7_CAGATC HXE64A2_CGATGT HXE64B8_ACTTGA HXL64A3_TTAGGC HXL64B9_GATCAG HXACA4_TGACCA HXACB10_TAGCTT HXECA5_ACAGTG HXECB11_GGCTAC HXLCA6_GCCAAT HXLCB12_CTTGTA HDA64A1_ATCACG HDA64B19_GTGAAA HDE64A4_TGACCA HDE64B20_GTGGCC HDL64A5_ACAGTG HDL64B22_CGTACG HDACA6_GCCAAT HDACB23_GAGTGG HDECA10_TAGCTT HDECB25_ACTGAT HDLCA12_CTTGTA HDLCB27_ATTCCT R19040024-HelainputD-Hl15R2hID_combined R19040024-HelainputD-Hl0h20JID_combined R19040024-HelainputD-Hl35R2hID_combined)

genomes=(genome_hg19_20kb_windows.bed genome_hg19_100kb_windows.bed genome_hg19_1Mb_windows.bed)


file=${samples[$SLURM_ARRAY_TASK_ID]}

getname="$(echo ${file} | sed 's/-/_/g')"
moreinfo="$(grep ${getname} ${mainPath}/project_repair_replication_all_samples.csv | sed 's/,/\t/g')"

if [ -d ${mainPath}/Input/$file ]; then
    samplePath=${mainPath}/Input/$file/pre_analysis
    minus="minus"
    plus="plus"
    mappedReads="$(grep -c '^' ${samplePath}/${file}_chr.bed)"
elif [ -d ${mainPath}/XRseq/$file ]; then
    samplePath=${mainPath}/XRseq/$file/pre_analysis
    minus="cutadapt_sorted_minus"
    plus="cutadapt_sorted_plus"
    mappedReads="$(grep -c '^' ${samplePath}/${file}_cutadapt_sorted_chr.bed)"
elif [ -d ${mainPath}/Damageseq/$file ]; then
    samplePath=${mainPath}/Damageseq/$file/pre_analysis
    minus="cutadapt_sorted_minus_dipyrimidines"
    plus="cutadapt_sorted_plus_dipyrimidines"
    minus_line="$(grep -c '^' ${samplePath}/${file}_cutadapt_sorted_minus_dipyrimidines.bed)"
    plus_line="$(grep -c '^' ${samplePath}/${file}_cutadapt_sorted_plus_dipyrimidines.bed)"
    mappedReads=`echo "$minus_line + $plus_line" | bc`
fi



for genome in ${genomes[@]}; do

bedtools intersect -a ${mainPath}/data/${genome} -b ${samplePath}/${file}_${minus}.bed -wa -c -F 0.5 > ${samplePath}/${file}_${genome%.*}_minus.txt

bedtools intersect -a ${mainPath}/data/${genome} -b ${samplePath}/${file}_${plus}.bed -wa -c -F 0.5 > ${samplePath}/${file}_${genome%.*}_plus.txt

${NGStoolkitPath}/addColumns.py -i ${samplePath}/${file}_${genome%.*}_minus.txt -o ${samplePath}/${file}_${genome%.*}_minus_full.txt -c "${moreinfo}" - "${mappedReads}" 

${NGStoolkitPath}/addColumns.py -i ${samplePath}/${file}_${genome%.*}_plus.txt -o ${samplePath}/${file}_${genome%.*}_plus_full.txt -c "${moreinfo}" + "${mappedReads}" 

cat ${samplePath}/${file}_${genome%.*}_plus_full.txt ${samplePath}/${file}_${genome%.*}_minus_full.txt > ${samplePath}/${file}_${genome%.*}_full.txt

python ${mainPath}/RPKM.py -i ${samplePath}/${file}_${genome%.*}_full.txt -c 4 -mr 0 -chse 2 3 -o ${samplePath}/${file}_${genome%.*}_full_rpkm.txt

done

echo "Finished!"          

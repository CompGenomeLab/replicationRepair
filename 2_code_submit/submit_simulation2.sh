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
#SBATCH --job-name=part2_sim
#SBATCH --account=users
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --qos=mid
#SBATCH --partition=mid
#SBATCH --time=23:00:00
#SBATCH --output=%j-slurm.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cemazgari@sabanciuniv.edu
#SBATCH --array=0-31

######13-16,29-35


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
#env
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


# sample list

samples=( NT4H1-CT5H2-HelaD1-5R2h2_combined XT4H1-CT5H2-HelaD3-5R2h2_combined R19029847-HD2TCD4-HelaD1-5R2h1_combined R19029847-HD2TCD4-HelaD3-5R2h1_combined HDA64A1_ATCACG HDA64B19_GTGAAA HDE64A4_TGACCA HDE64B20_GTGGCC HDL64A5_ACAGTG HDL64B22_CGTACG HDACA6_GCCAAT HDACB23_GAGTGG HDECA10_TAGCTT HDECB25_ACTGAT HDLCA12_CTTGTA HDLCB27_ATTCCT HXA64A1_ATCACG HXA64B7_CAGATC HXE64A2_CGATGT HXE64B8_ACTTGA HXL64A3_TTAGGC HXL64B9_GATCAG HXACA4_TGACCA HXACB10_TAGCTT HXECA5_ACAGTG HXECB11_GGCTAC HXLCA6_GCCAAT HXLCB12_CTTGTA R19026421-2019XR1-Hela15X2_combined_R1 R19026421-2019XR1-Hela35X3_combined_R1 R19033030-2019XR3-Hela15X7_combined_R1 R19033030-2019XR3-Hela35X8_combined_R1 )


file=${samples[$SLURM_ARRAY_TASK_ID]}

echo This is task $SLURM_ARRAY_TASK_ID

mainPath="/cta/users/cemazgari/repairReplication"
countPath="${mainPath}/git/replicationRepair/4_output/data"
NGStoolkitPath="/cta/users/cemazgari/NGStoolkit/bin"
getname="$(echo $file | sed 's/-/_/g')"
sample="$(grep ${getname} ${mainPath}/project_repair_replication_all_samples.csv | sed 's/,/\t/g' | awk '{print $1}')"
counts="$(grep ${sample} ${countPath}/read_counts.csv | sed 's/,/\t/g' | awk '{print $2}')"

echo $file, $sample, $counts

if [[ $sample == *"XR"*"120"* ]]; then 

${mainPath}/simulation/filter_syn_fasta -in ${mainPath}/simulation/inputfasta/${file}_cutadapt_sorted_chr.fa -length 26 -numReads ${counts} -sim ${mainPath}/simulation/XR_sim_f3.fq -out ${mainPath}/simulation/rawdata/${file}_sim.fastq

elif [[ $sample == *"XR"*"12_"* ]]; then

${mainPath}/simulation/filter_syn_fasta -in ${mainPath}/simulation/inputfasta/${file}_cutadapt_sorted_chr.fa -length 26 -numReads ${counts} -sim ${mainPath}/simulation/XR_sim_f1.fq -out ${mainPath}/simulation/rawdata/${file}_sim.fastq

elif [[ $sample == *"DS"*"120"* ]]; then 

${mainPath}/simulation/filter_syn_fasta -in ${mainPath}/simulation/inputfasta/${file}_cutadapt_sorted_dipyrimidines.fa -length 10 -numReads ${counts} -sim ${mainPath}/simulation/DS_sim_f3.fq -out ${mainPath}/simulation/rawdata/${file}_sim.fastq

elif [[ $sample == *"DS"*"12_"* ]]; then

${mainPath}/simulation/filter_syn_fasta -in ${mainPath}/simulation/inputfasta/${file}_cutadapt_sorted_dipyrimidines.fa -length 10 -numReads ${counts} -sim ${mainPath}/simulation/DS_sim_f1.fq -out ${mainPath}/simulation/rawdata/${file}_sim.fastq

fi

awk 'NR%4==1{a=substr($0,2);}NR%4==2{print ">"a"\n"$0}' ${mainPath}/simulation/rawdata/${file}_sim.fastq > ${mainPath}/simulation/rawdata/${file}_sim.fasta

${NGStoolkitPath}/fa2kmerAbundanceTable.py -i ${mainPath}/simulation/rawdata/${file}_sim.fasta -k 2 -o ${mainPath}/simulation/rawdata/${file}_sim_dinucleotideTable.txt
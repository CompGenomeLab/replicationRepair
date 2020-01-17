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
#SBATCH --job-name=Hela_input
#SBATCH --account=users
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --qos=long
#SBATCH --partition=long
#SBATCH --time=7-00
#SBATCH --output=%j-slurm.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cemazgari@sabanciuniv.edu
#SBATCH --array=0-2



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
module load trimmomatic-0.36-gcc-8.2.0-cw2gsdk
module load salmon-0.9.1-gcc-7.2.0-mlvecrk

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

rawdataPath="/cta/groups/adebali/data/repair/191005"
mainPath="/cta/users/cemazgari/repairReplication"
genomePath="/cta/groups/adebali/data/reference_genomes/human/gencode/19"

samples=(R19040024-HelainputD-Hl0h20JID_combined R19040024-HelainputD-Hl15R2hID_combined R19040024-HelainputD-Hl35R2hID_combined)

file=${samples[$SLURM_ARRAY_TASK_ID]}

echo $file

mkdir -p ${mainPath}/Input/${file}

bowtie2 -X 1000 -x ${genomePath}/Bowtie2/genome -1 ${rawdataPath}/${file}_R1.fastq.gz -2 ${rawdataPath}/${file}_R2.fastq.gz -S ${mainPath}/Input/${file}/${file}.sam

samtools view -q 20 -b -o ${mainPath}/Input/${file}/${file}.bam ${mainPath}/Input/${file}/${file}.sam # samtools: sam to bam

bedtools bamtobed -i ${mainPath}/Input/${file}/${file}.bam > ${mainPath}/Input/${file}/${file}.bed # bedtools: bam to bed

grep chr ${mainPath}/Input/${file}/${file}.bed > ${mainPath}/Input/${file}/${file}_chr.bed

awk '$6=="-"{print $0}' ${mainPath}/Input/${file}/${file}_chr.bed | sort -k1,1 -k2,2n -k3,3n > ${mainPath}/Input/${file}/${file}_minus.bed

awk '$6=="+"{print $0}' ${mainPath}/Input/${file}/${file}_chr.bed | sort -k1,1 -k2,2n -k3,3n > ${mainPath}/Input/${file}/${file}_plus.bed

echo "done!"




 





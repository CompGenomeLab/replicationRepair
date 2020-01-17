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
#SBATCH --job-name=simulation_part1
#SBATCH --account=users
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --qos=long
#SBATCH --partition=long
#SBATCH --time=6-23
#SBATCH --output=%j-slurm.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cemazgari@sabanciuniv.edu



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

#main dir for simulation
outPath="/cta/users/cemazgari/repairReplication/simulation"
input="/cta/groups/adebali/data/reference_genomes/human/gencode/19/genome.fa"
toolPath="/cta/users/cemazgari/repairReplication/simulation/art_bin_GreatSmokyMountains"

${toolPath}/art_illumina -ss HS25 -i ${input} -l 10 -f 10 -o ${outPath}/human_sim_f10 -rs 1002

#!/bin/bash

source ../1_code/source_dir.sh

# Download genome

wget ${genomeDownload}/release_19/GRCh37.p13.genome.fa.gz -P ${genomePath}/ 

gunzip ${genomePath}/GRCh37.p13.genome.fa.gz

# Preparing indexes

mv ${genomePath}/GRCh37.p13.genome.fa ${genomePath}/genome.fa 

bowtie2-build ${genomePath}/genome.fa  ${genomePath}/Bowtie2/genome

samtools faidx ${genomePath}/genome.fa 

#



#!/bin/bash

	genomeDownload="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
	genomePath="/cta/groups/adebali/data/reference_genomes/human/gencode"

	# Human Genome Release 29

		wget ${genomeDownload}/release_29/GRCh38.p12.genome.fa.gz -P ${genomePath}/29/ # Must be done on login node

		gunzip ${genomePath}/29/GRCh38.p12.genome.fa.gz


	# Human Genome Release 19
		
		wget ${genomeDownload}/release_19/GRCh37.p13.genome.fa.gz -P ${genomePath}/19/ # Must be done on login node

		gunzip ${genomePath}/19/GRCh37.p13.genome.fa.gz

    #


# Preparing indexes

	# Human Genome Release 29

		mv ${genomePath}/29/GRCh38.p12.genome.fa ${genomePath}/29/genome.fa 
		
		bowtie2-build ${genomePath}/29/genome.fa  ${genomePath}/29/Bowtie2/genome

		samtools faidx ${genomePath}/29/genome.fa 

	# Human Genome Release 19

		mv ${genomePath}/19/GRCh37.p13.genome.fa ${genomePath}/19/genome.fa 

		bowtie2-build ${genomePath}/19/genome.fa  ${genomePath}/19/Bowtie2/genome

		samtools faidx ${genomePath}/19/genome.fa 

	#


# 

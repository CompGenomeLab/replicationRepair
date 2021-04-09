#!/bin/bash

# Please provide the full directory of the repository
repositoryPath=""  # do not leave "/" at the end of the path

mainPath="${repositoryPath}/repairReplication" 
rawdataPath="${mainPath}/0_data/raw"
simRawPath="${mainPath}/0_data/simulation"
codePath="${mainPath}/1_code"
genomePath="${mainPath}/0_data/hg19"
genomeDownload="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
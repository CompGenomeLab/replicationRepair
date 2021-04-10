# Divergent effects of DNA replication on the repair of UV-induced DNA damage. 

This repository contains the codes that are used for the analyses of the paper "_Divergent effects of DNA replication on the repair of UV-induced DNA damage._"

## Authors

_Yanchao Huang, Cem Azgari, Yi-Ying Chiou, Laura A. Lindsey-Boltz, Aziz Sancar, Jinchuan Hu, Ogun Adebali_

## Description 

UV-induced damage can cause mutations during DNA replication. The crosstalk between replication and repair may contribute to UV-related mutagenesis. By integrating genome-wide damage and repair maps along with replication maps, we investigated the effects of DNA replication on nucleotide excision repair. Early replication domains are repaired faster due to open chromatin; thus they harbor fewer mutations. Ongoing replication can exert an additional impact of promoting local repair by relaxing surrounding chromatin. Repair levels also show a strand asymmetry favoring leading strands, presumably because leading strand synthesis is more active to recover double-strand after encountering a lesion. This biased repair coincides with the replicative mutation asymmetry in melanoma, indicating a role of exogenous damage and repair in replication-associated mutation asymmetry.

### Required Programs

- cutadapt
- bowtie2
- bedtools
- bamtools
- samtools
- [ART](https://pubmed.ncbi.nlm.nih.gov/22199392/): a next-generation sequencing read simulator

### Usage

* Clone the repository and navigate into the directory: 

```
git clone https://github.com/CompGenomeLab/replicationRepair.git
    
cd replicationRepair
```

* Raw data that is analyzed in this repository can be found with the accession code PRJNA608124 and can be downloaded by `sra-toolkit`. Any raw data must be stored at `0_data/raw/` directory.

* Initially, genome files must be downloaded to run the pipelines, which can be done by using `mypre_repairReplication.sh` pipeline:

```
./2_pipeline/mypre_repairReplication.sh
```


| Warning: please don't forget to give execution permission to the pipeline files. |
| --- |


* The output genome files will be stored at `/0_data/hg19`. To process the XR-seq and Damage-seq samples, the full directory of the repository must be given to `source_dir.sh` script (located at `1_code/`) and can be run only providing the name of the sample without the extension:

```
./2_pipeline/myXRseq.sh ${XRseqSample}

./2_pipeline/mydamageseq.sh ${DamageseqSample}
```

| Warning: if the layout is paired, you will have 2 samples. In that case, if sample names are `sample_(R1/1).fastq.gz` and `sample_(R2/2).fastq.gz`, you can provide the name as `sample`. The pipeline will notice that there are 2 files in the directory, thus will set the layout as paired. |
| --- |

* All the details about the XR-seq and Damage-seq pipelines can be found at `2_pipelines/` directory.

* After the pipelines are processed they should be reorganized with `1_final_report.R` script and each script for figures, which are tagged with the corresponding figure number, can be found in `1_code/` directory.

* We have also provided the processed files to generate the figures and the figures themselves, which can be found at `3_output/processed_files/` and `3_output/plots`, respectively.
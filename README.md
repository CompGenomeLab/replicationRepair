# Divergent effects of DNA replication on the repair of UV-induced DNA damage. 

This repository contains the codes that are used for the analyses of the paper "_Divergent effects of DNA replication on the repair of UV-induced DNA damage._"

## Authors

_Yanchao Huang, Cem Azgari, Yi-Ying Chiou, Laura A. Lindsey-Boltz, Aziz Sancar, Jinchuan Hu, Ogun Adebali_

## Description 

UV-induced damage can cause mutations during DNA replication. The crosstalk between replication and repair may contribute to UV-related mutagenesis. By integrating genome-wide damage and repair maps along with replication maps, we investigated the effects of DNA replication on nucleotide excision repair. Early replication domains are repaired faster due to open chromatin; thus they harbor fewer mutations. Ongoing replication can exert an additional impact of promoting local repair by relaxing surrounding chromatin. Repair levels also show a strand asymmetry favoring leading strands, presumably because leading strand synthesis is more active to recover double-strand after encountering a lesion. This biased repair coincides with the replicative mutation asymmetry in melanoma, indicating a role of exogenous damage and repair in replication-associated mutation asymmetry.

## Installation

- This workflow is prepared using 
[Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management 
system and [conda](https://docs.conda.io/en/latest/)

- To run the workflow, you should have conda installed for environment 
management. All the other packages including Snakemake and their dependencies 
can be obtained automatically through environments prepared for each step of 
the workflow. You can follow the installation steps from 
[the link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html).

- Initially, you should clone the repository and navigate into the directory: 

    ```
    git clone https://github.com/CompGenomeLab/replicationRepair.git
    
    cd replicationRepair
    ```

- Next, you should create a conda environment with the defined packages. 
Install [mamba](https://mamba.readthedocs.io/en/latest/) 
and create the environment using mamba:

    ```
    conda install -c conda-forge mamba

    mamba create -c bioconda -c conda-forge -c r -n repair snakemake=6.3.0 python=3.8 rust=1.50 sra-tools=2.11.0

    conda activate repair
    ```

<br>

## Directory Structure

This workflow is prepared according to the 
[structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) 
recommended by Snakemake: 

- `config/`: contains the configuration files.

- `logs/`: contains the log files of each step. 
This folder will automatically appear when you run the workflow.

- `reports/`: contains the report files, which can be produced 
after the workflow is over. 

- `resources/`: contains `samples/` where the raw XR-seq and Damage-seq data 
are stored and `ref_genomes/` where the reference genome files are stored. 

- `results/`: contains the generated files and figures. *** DETAYLI ANLATIM

- `workflow/`: contains `envs/` where the environments are stored, 
`rules/` where the Snakemake rules are stored, 
`scripts/` where the scripts used inside the rules are stored, and
`snakefiles/` where rule flow of each pipeline of the project can be found.
<br>

## Usage

Before *** XR damage çalışacak, hg19 referansı eklenmeli

After adjusting the configuration file, you can run the workflow 
from `replicationRepair` directory:

    ```
    snakemake --cores 64 --use-conda --keep-going
    ```

| Note: To run the workflow on [Slurm Workload Manager](https://slurm.schedmd.com/srun.html) as set of jobs, `--profile` flag must be provided with proper slurm configuration file (`config/slurm`). |
| --- |









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

* To process melanoma mutations data, `mutation_analysis.sh` pipeline must run.

```
./2_pipeline/mutation_analysis.sh
```

* Note: Before running the pipeline, the name of the mutation file and the data file, which mutations will be mapped on, must be provided (with full path).  

* After the pipelines are processed they should be reorganized with `1_final_report.R` script and each script for figures, which are tagged with the corresponding figure number, can be found in `1_code/` directory.

* We have also provided the processed files to generate the figures and the figures themselves, which can be found at `3_output/processed_files/` and `3_output/plots`, respectively.
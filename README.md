# The effects of replication domains on the genome-wide UV-induced DNA damage and repair 

This repository contains the codes that are used for the analyses of the paper "_The effects of replication domains on the genome-wide UV-induced DNA damage and repair_"

## Authors

_Yanchao Huang, Cem Azgari, Mengdie Yin, Yi-Ying Chiou, Laura A. Lindsey-Boltz, Aziz Sancar, Jinchuan Hu, Ogun Adebali_

## Description 

Nucleotide excision repair is the primary repair mechanism that removes UV-induced DNA lesions in placentals. If the UV-induced lesions are left unrepaired they might turn into mutations during DNA replication. Although the mutagenesis of pyrimidine dimers is reasonably well understood, the direct effects of replication fork progress on nucleotide excision repair are yet to be clarified. Here, we applied Damage-seq and XR-seq techniques and generated replication maps in synchronized UV-treated HeLa cells. The results suggested that ongoing replication stimulates local repair in both early and late replication domains by relaxing surrounding chromatin. On the other hand, it was unveiled that lesions on lagging strand templates were repaired slower in late replication domains, which was probably due to the imbalanced sequence context. The asymmetric relative repair was in line with the strand bias of melanoma mutations, suggesting a role of exogenous damage, repair, and replication in the mutational strand asymmetry.

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

    mamba create -c bioconda -c conda-forge -c r -n repair snakemake=6.3.0

    conda activate repair
    ```

### Genome Download

- Genome fasta file of hg19 can be downloaded and unzipped with the commands below:

    ```
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz

    gunzip GRCh37.p13.genome.fa.gz
    ```

- Then the fasta file should be named as `genome_hg19.fa` and 
moved to `resources/ref_genomes/hg19/` directory located in 
`replicationRepair`.

### XR-seq and Damage-seq Pipelines

- XR-seq and Damage-seq pipelines should be cloned separately 
(outside of `replicationRepair` directory) from the 
[github link](https://github.com/CompGenomeLab/xr-ds-seq-snakemake).

- To produce consistent results, `genome_hg19.fa` should be copied to
`xr-ds-seq-snakemake/resources/ref_genomes/hg19/` 
(`hg19/` directory should be created before copying the genome file).

- Lastly, `config_DS.yaml` and `config_XR.yaml` files in `xr-ds-seq-snakemake/config/` 
should be replaced by their counterparts in `replicationRepair/config/`.

### Retrieve Melanoma Simple Somatic Mutations 

- Simple somatic mutations of melanoma are publicly available in 
[ICGC Data Portal](https://dcc.icgc.org/releases/release_28/Projects/MELA-AU).
You should download `simple_somatic_mutation.open.MELA-AU.tsv.gz` file, 
rename it as `melanoma.tsv.gz`, and 
move it to the `replicationRepair/resources/samples/mutation` directory.

<br>

## Directory Structure

This workflow is prepared according to the 
[structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) 
recommended by Snakemake: 

- `config/`: contains the configuration files.

- `logs/`: contains the log files of each step. 
This folder will automatically appear when you run the workflow.

- `report/`: contains the description files of figures,
which will be used in reports.

- `resources/`: contains `samples/` where the raw XR-seq and Damage-seq data 
are stored and `ref_genomes/` where the reference genome files are stored. 

- `results/`: contains the generated files and figures.

- `workflow/`: contains `envs/` where the environments are stored, 
`rules/` where the Snakemake rules are stored, and
`scripts/` where the scripts used inside the rules are stored.
<br>

## Usage

### XR-seq and Damage-seq

- Initially the XR-seq and Damage-seq pipeline should be run
with the provided config files. 

- After the pipeline is competed, produced bed files 
should be moved to the appropriate directories.

    - For XR-seq samples:
    ```
    cp {path_to_dir}/xr-ds-seq-snakemake/results/XR/*/*_*us.bed {path_to_dir}/replicationRepair/resources/samples/XR/
    ```

    - For Damage-seq samples:
    ```
    cp {path_to_dir}/xr-ds-seq-snakemake/results/DS/*/*_ds_dipyrimidines_*us.bed {path_to_dir}/replicationRepair/resources/samples/DS/
    ```

    - For simulated samples:
    ```
    cp {path_to_dir}/xr-ds-seq-snakemake/results/*/*/*_sim.bed {path_to_dir}/replicationRepair/resources/samples/sim/
    ```

### Run Further Analyses

You can run the workflow from `replicationRepair` directory:

    snakemake --cores 64 --use-conda --keep-going

| Note: To run the workflow on [Slurm Workload Manager](https://slurm.schedmd.com/srun.html) as set of jobs, `--profile` flag must be provided with proper slurm configuration file (`config/slurm`). |
| --- |

<br>

To generate detailed HTML report files, 
the code below should be run after workflow:

```
snakemake --report report.zip
```
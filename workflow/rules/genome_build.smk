
rule genome_build:
    input:
        reference="resources/ref_genomes/hg19/genome_hg19.fa",
    output:
        multiext(
        "resources/ref_genomes/hg19/Bowtie2/genome_hg19",
        ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    log: 
        "logs/rule/analysis/genome_build.log"
    benchmark:
        "logs/rule/analysis/genome_build.benchmark.txt",
    params:
        extra=""  
    threads: 
        4
    wrapper:
        "0.69.0/bio/bowtie2/build"
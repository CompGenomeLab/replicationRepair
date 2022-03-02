
rule genome_indexing:
    input:
        "resources/ref_genomes/hg19/genome_hg19.fa",
    output:
        "resources/ref_genomes/hg19/genome_hg19.fa.fai",
    benchmark:
        "logs/rule/analysis/indexing.benchmark.txt",
    wrapper:
        "0.69.0/bio/samtools/faidx"
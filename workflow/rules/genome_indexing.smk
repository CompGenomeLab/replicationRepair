
rule genome_indexing:
    input:
        "resources/ref_genomes/hg19/genome_hg19.fa",
    output:
        "resources/ref_genomes/hg19/genome_hg19.fa.fai",
    benchmark:
        "logs/rule/analysis/indexing.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:
        """
        (echo "`date -R`: Creating fai file..." &&
        samtools faidx {input} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """
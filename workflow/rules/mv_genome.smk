rule mv_genome:
    input:
        "resources/ref_genomes/hg19/genome_hg19_50kb.bed",
    output:
        "results/regions/genome_hg19_50kb.bed",
    log:
        "logs/rule/analysis/regions/mv_genome.log",
    benchmark:
        "logs/rule/analysis/regions/mv_genome.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Make bed file stranded..." &&
        awk '{{print $1"\\t"$2"\\t"$3"\\t"".""\\t"".""\\t""."}}' {input} \
        > {output} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

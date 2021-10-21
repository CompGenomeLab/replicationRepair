rule intergenic_mutation:
    input:
        plus="results/mutation/{samples}/{samples}_target_mut_plus.tsv",
        minus="results/mutation/{samples}/{samples}_target_mut_minus.tsv",
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes.bed",
    output:
        plus_intergenic="results/mutation/{samples}/{samples}_intergenic_target_mut_plus.tsv",
        minus_intergenic="results/mutation/{samples}/{samples}_intergenic_target_mut_minus.tsv",
    log:
        "logs/{samples}/{samples}_intergenic_mutation.log",
    benchmark:
        "logs/{samples}/{samples}_intergenic_mutation.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Getting intergenic (plus strand)..." &&
        bedtools intersect \
        -a {input.plus} \
        -b {input.genes} \
        -v -f 0.5 > {output.plus_intergenic} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Getting intergenic (minus strand)..." &&
        bedtools intersect \
        -a {input.minus} \
        -b {input.genes} \
        -v -f 0.5 > {output.minus_intergenic} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
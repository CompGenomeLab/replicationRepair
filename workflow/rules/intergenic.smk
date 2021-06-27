rule intergenic_xr:
    input:
        plus="resources/samples/XR/{samples}_{build}_sorted_plus.bed",
        minus="resources/samples/XR/{samples}_{build}_sorted_minus.bed",
        genes="results/regions/hg19_ucsc_genes_knownCanonical_stranded.bed",
    output:
        plus_intergenic=temp("results/XR/{samples}/{samples}_{build}_intergenic_sorted_plus.bed"),
        minus_intergenic=temp("results/XR/{samples}/{samples}_{build}_intergenic_sorted_minus.bed"),
    log:
        "logs/{samples}/{samples}_{build}_intergenic_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intergenic_xr.benchmark.txt",
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
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Getting intergenic (minus strand)..." &&
        bedtools intersect \
        -a {input.minus} \
        -b {input.genes} \
        -v -f 0.5 > {output.minus_intergenic} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intergenic_ds:
    input:
        plus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
        genes="results/regions/hg19_ucsc_genes_knownCanonical_stranded.bed",
    output:
        plus_intergenic=temp("results/DS/{samples}/{samples}_{build}_intergenic_sorted_ds_dipyrimidines_plus.bed"),
        minus_intergenic=temp("results/DS/{samples}/{samples}_{build}_intergenic_sorted_ds_dipyrimidines_minus.bed"), 
    log:
        "logs/{samples}/{samples}_{build}_intergenic_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intergenic_ds.benchmark.txt",
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
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Getting intergenic (minus strand)..." &&
        bedtools intersect \
        -a {input.minus} \
        -b {input.genes} \
        -v -f 0.5 > {output.minus_intergenic} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intergenic_sim:
    input:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus.bed",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus.bed",
        genes="results/regions/hg19_ucsc_genes_knownCanonical_stranded.bed",
    output:
        plus_intergenic=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_plus_intergenic.bed"),
        minus_intergenic=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_minus_intergenic.bed"),
    log:
        "logs/{samples}/{samples}_{build}_intergenic_{method}_sim.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intergenic_{method}_sim.benchmark.txt",
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
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Getting intergenic (minus strand)..." &&
        bedtools intersect \
        -a {input.minus} \
        -b {input.genes} \
        -v -f 0.5 > {output.minus_intergenic} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

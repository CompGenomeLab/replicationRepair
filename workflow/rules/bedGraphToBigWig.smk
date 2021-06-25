
rule bedGraphToBigWig_input:
    input:
        bdg="results/input/{samples}/{samples}_{build}_sorted_{strand}.bdg",
        index="resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        report("results/input/{samples}/{samples}_{build}_sorted_{strand}.bw", 
                category="BigWig"),
    log:
        "logs/{samples}/{samples}_{build}_bedGraphToBigWig_{strand}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bedGraphToBigWig_{strand}.benchmark.txt",
    conda:
        "../envs/bedGraphToBigWig.yaml"
    shell:  
        """
        (echo "`date -R`: Converting bedGraph to bigWig..." &&
        bedGraphToBigWig {input.bdg} {input.index} {output} &&
        echo "`date -R`: Success! Conversion is done." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """
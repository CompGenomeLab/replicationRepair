rule bedGraphToBigWig_edu:
    input:
        bdg="results/edu/{samples}/{samples}_{build}_sorted_{strand}.bdg",
        index="resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        "results/edu/{samples}/{samples}_{build}_sorted_{strand}.bw",
    log:
        "logs/{samples}/{samples}_{build}_bedGraphToBigWig_{strand}_edu.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bedGraphToBigWig_{strand}_edu.benchmark.txt",
    conda:
        "../envs/bedGraphToBigWig.yaml"
    shell:  
        """
        (echo "`date -R`: Converting bedGraph to bigWig..." &&
        bedGraphToBigWig {input.bdg} {input.index} {output} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
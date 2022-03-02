rule bedGraphToBigWig_edu:
    input:
        bdg="results/edu/{samples}/{samples}_hg19_sorted_{strand}.bdg",
        index="resources/ref_genomes/hg19/genome_hg19.fa.fai",
    output:
        "results/edu/{samples}/{samples}_hg19_sorted_{strand}.bw",
    log:
        "logs/rule/analysis/{samples}/{samples}_hg19_bedGraphToBigWig_{strand}_edu.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_hg19_bedGraphToBigWig_{strand}_edu.benchmark.txt",
    conda:
        "../envs/bedGraphToBigWig.yaml"
    shell:  
        """
        (echo "`date -R`: Converting bedGraph to bigWig..." &&
        bedGraphToBigWig {input.bdg} {input.index} {output} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
rule produceInitiationZones:
    input:  
        okseq="results/okseq/{samples}/{samples}_se_{build}_sorted.bam",  
        genome="resources/ref_genomes/{build}/genome_{build}.txt",           
    output: 
        "results/okseq/{samples}/{samples}_{build}_HMMsegments_IZ.bed",  
    params:
        prefix="results/okseq/{samples}/{samples}_{build}",
    log:
        "logs/{samples}_{build}_produceInitiationZones.log",
    benchmark:
        "logs/{samples}_{build}_produceInitiationZones.benchmark.txt",
    conda:
        "../envs/produceInitiationZones.yaml",
    shell:
        """
        Rscript workflow/scripts/okseqhmm.R \
        -i {input.okseq} \
        --genome {input.genome} \
        -p {params.prefix} \
        -o {output} &> {log}
        """
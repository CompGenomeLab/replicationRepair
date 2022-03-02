rule produceInitiationZones:
    input:  
        okseq="results/okseq/{samples}/{samples}_se_hg19_sorted.bam",  
        genome="resources/ref_genomes/hg19/genome_hg19.txt",           
    output: 
        "results/okseq/{samples}/{samples}_hg19_HMMsegments_IZ.bed",  
    params:
        prefix="results/okseq/{samples}/{samples}_hg19",
    log:
        "logs/rule/analysis/{samples}_hg19_produceInitiationZones.log",
    benchmark:
        "logs/rule/analysis/{samples}_hg19_produceInitiationZones.benchmark.txt",
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
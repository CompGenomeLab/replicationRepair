
rule cutadapt_se_okseq:
    input:
        "resources/samples/okseq/{samples}.fastq.gz",
    output:
        fastq=temp("results/okseq/{samples}/{samples}_cutadapt.fastq.gz"),
        qc=report("results/okseq/{samples}/{samples}_cutadapt.qc.txt", category="QC"),   
    params:
        adapters="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",  
        extra="",  
    threads: 32
    log:
        "logs/rule/cutadapt_se_okseq/{samples}.log",
    benchmark:
        "logs/rule/cutadapt_se_okseq/{samples}.benchmark.txt",
    conda:
        "../envs/cutadapt.yaml"
    shell:
        """
        (echo "`date -R`: Trimming adapters..." &&
        cutadapt \
        -j {threads} \
        {params.adapters} \
        {params.extra} \
        -o {output.fastq} {input} \
        > {output.qc} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """


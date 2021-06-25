
rule cutadapt_se_okseq:
    input:
        "resources/samples/okseq/{samples}.fastq.gz",
    output:
        fastq=temp("results/okseq/{samples}/{samples}_cutadapt.fastq.gz"),
        qc=report("results/okseq/{samples}/{samples}_cutadapt.qc.txt", category="QC"),   
    params:
        adapters="-a ACACTCTTTCCCTACACGACGCTCTTCC",  
        extra="",  
    log:
        "logs/{samples}/{samples}_cutadapt.log",
    benchmark:
        "logs/{samples}/{samples}_cutadapt.benchmark.txt",
    wrapper:
        "0.69.0/bio/cutadapt/se"

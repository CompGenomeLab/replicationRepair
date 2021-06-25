
rule fastqc_se_input:
    input:
        "resources/samples/input/{samples}.fastq.gz",
    output:
        html=report("results/input/{samples}/{samples}.html", category="QC"),
        zip="results/input/{samples}/{samples}_fastqc.zip",
    params: ""
    log:
        "logs/{samples}/{samples}_fastqc.log",
    benchmark:
        "logs/{samples}/{samples}_fastqc.benchmark.txt",
    threads: 1
    wrapper:
        "0.69.0/bio/fastqc"

rule fastqc_pe_input:
    input:
        "resources/samples/input/{samples}_{ext}.fastq.gz", 
    output:
        html=report("results/input/{samples}/{samples}_{ext}.html", category="QC"), 
        zip="results/input/{samples}/{samples}_{ext}_fastqc.zip", 
    params: ""
    log:
        "logs/{samples}/{samples}_fastqc_{ext}.log", 
    benchmark:
        "logs/{samples}/{samples}_fastqc_{ext}.benchmark.txt",
    threads: 1
    wrapper:
        "0.69.0/bio/fastqc"

rule fastqc_se_okseq:
    input:
        "resources/samples/okseq/{samples}.fastq.gz",
    output:
        html=report("results/okseq/{samples}/{samples}.html", category="QC"),
        zip="results/okseq/{samples}/{samples}_fastqc.zip",
    params: ""
    log:
        "logs/{samples}/{samples}_fastqc.log",
    benchmark:
        "logs/{samples}/{samples}_fastqc.benchmark.txt",
    threads: 1
    wrapper:
        "0.69.0/bio/fastqc"

rule fastqc_pe_okseq:
    input:
        "resources/samples/okseq/{samples}_{ext}.fastq.gz", 
    output:
        html=report("results/okseq/{samples}/{samples}_{ext}.html", category="QC"), 
        zip="results/okseq/{samples}/{samples}_{ext}_fastqc.zip", 
    params: ""
    log:
        "logs/{samples}/{samples}_fastqc_{ext}.log", 
    benchmark:
        "logs/{samples}/{samples}_fastqc_{ext}.benchmark.txt",
    threads: 1
    wrapper:
        "0.69.0/bio/fastqc"
rule fastqc_se_okseq:
    input:
        "resources/samples/okseq/{samples}.fastq.gz",
    output:
        html=report("results/okseq/{samples}/{samples}.html", category="QC"),
        zip="results/okseq/{samples}/{samples}_fastqc.zip",
    params: ""
    log:
        "logs/rule/analysis/{samples}/{samples}_fastqc_se_okseq.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_fastqc_se_okseq.benchmark.txt",
    threads: 16
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
        "logs/rule/analysis/{samples}/{samples}_fastqc_pe_okseq_{ext}.log", 
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_fastqc_pe_okseq_{ext}.benchmark.txt",
    threads: 16
    wrapper:
        "0.69.0/bio/fastqc"

rule fastqc_se_edu:
    input:
        "resources/samples/edu/{samples}.fastq.gz",
    output:
        html=report("results/edu/{samples}/{samples}.html", category="QC"),
        zip="results/edu/{samples}/{samples}_fastqc.zip",
    params: ""
    log:
        "logs/rule/analysis/{samples}/{samples}_fastqc_se_edu.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_fastqc_se_edu.benchmark.txt",
    threads: 16
    wrapper:
        "0.69.0/bio/fastqc"

rule fastqc_pe_edu:
    input:
        "resources/samples/edu/{samples}_{ext}.fastq.gz", 
    output:
        html=report("results/edu/{samples}/{samples}_{ext}.html", category="QC"), 
        zip="results/edu/{samples}/{samples}_{ext}_fastqc.zip", 
    params: ""
    log:
        "logs/rule/analysis/{samples}/{samples}_fastqc_pe_edu_{ext}.log", 
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_fastqc_pe_edu_{ext}.benchmark.txt",
    threads: 16
    wrapper:
        "0.69.0/bio/fastqc"
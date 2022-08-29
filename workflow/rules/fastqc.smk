rule fastqc_se_okseq:
    input:
        "resources/samples/okseq/{samples}.fastq.gz",
    output:
        html=report("results/okseq/{samples}/{samples}_fastqc.html", category="QC"),
        zip="results/okseq/{samples}/{samples}_fastqc.zip",
    params: 
        extra="",
        tmpdir="results/okseq/{samples}",
    log:
        "logs/rule/fastqc_se_okseq/{samples}.log",
    benchmark:
        "logs/rule/fastqc_se_okseq/{samples}.benchmark.txt",
    threads: 16
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        (echo "`date -R`: FastQC..." &&
        fastqc \
        {params.extra} \
        -t {threads} \
        --outdir {params.tmpdir} \
        {input} && 
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1 
        """

rule fastqc_se_edu:
    input:
        "resources/samples/edu/{samples}.fastq.gz",
    output:
        html=report("results/edu/{samples}/{samples}_fastqc.html", category="QC"),
        zip="results/edu/{samples}/{samples}_fastqc.zip",
    params: 
        extra="",
        tmpdir="results/edu/{samples}",
    log:
        "logs/rule/fastqc_se_edu/{samples}.log",
    benchmark:
        "logs/rule/fastqc_se_edu/{samples}.benchmark.txt",
    threads: 16
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        (echo "`date -R`: FastQC..." &&
        fastqc \
        {params.extra} \
        -t {threads} \
        --outdir {params.tmpdir} \
        {input} && 
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1 
        """

rule fastqc_pe_edu:
    input:
        "resources/samples/edu/{samples}_{ext}.fastq.gz", 
    output:
        html=report("results/edu/{samples}/{samples}_{ext}_fastqc.html", category="QC"), 
        zip="results/edu/{samples}/{samples}_{ext}_fastqc.zip", 
    params: 
        extra="",
        tmpdir="results/edu/{samples}_{ext}",
        name="{samples}_{ext}",
    log:
        "logs/rule/fastqc_pe/{samples}_{ext}.log",
    benchmark:
        "logs/rule/fastqc_pe/{samples}_{ext}.benchmark.txt",
    threads: 16
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        (echo "`date -R`: FastQC..." &&
        mkdir -p {params.tmpdir} &&
        fastqc \
        {params.extra} \
        -t {threads} \
        --outdir {params.tmpdir} \
        {input} &&
        mv {params.tmpdir}/{params.name}_fastqc.html {output.html} && 
        mv {params.tmpdir}/{params.name}_fastqc.zip {output.zip} &&
        rmdir {params.tmpdir} && 
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1 
        """


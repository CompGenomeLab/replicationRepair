rule unzipTSS:
    input:
        "resources/ref_genomes/hg19/hg19_ucsc_genes_knownCanonical_{tss_tes}_windows_201_100.bed.gz",
    output:
        "resources/ref_genomes/hg19/hg19_ucsc_genes_knownCanonical_{tss_tes}_windows_201_100.bed",
    log:
        "logs/rule/analysis/unzipTSS_{tss_tes}.log",
    benchmark:
        "logs/rule/analysis/unzipTSS_{tss_tes}.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Unzip {input}..." &&
        gunzip {input} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
rule genome_download:
    output:
        "resources/ref_genomes/hg19/genome_hg19.fa"
    params:
        link=config["genome"]["link"],
    log:
        "logs/rule/genome_download/hg19.log",
    benchmark:
        "logs/rule/genome_download/hg19.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Downloading hg19 genome..." &&
        wget {params.link} -O {output}.gz &&
        gunzip {output}.gz &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
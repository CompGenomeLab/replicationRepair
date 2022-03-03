
rule genome_build:
    input:
        reference="resources/ref_genomes/hg19/genome_hg19.fa",
    output:
        multiext(
        "resources/ref_genomes/hg19/Bowtie2/genome_hg19",
        ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    log: 
        "resources/ref_genomes/hg19/log/bowtie2_build.log"
    benchmark:
        "resources/ref_genomes/hg19/log/bowtie2_build.benchmark.txt",
    params:
        extra="", 
        name="resources/ref_genomes/hg19/Bowtie2/genome_hg19",
    conda:
        "../envs/align.yaml" 
    threads: 
        16
    shell: 
        """
        (echo "`date -R`: Building indexes..." &&
        bowtie2-build --threads {threads} \
        {params.extra} \
        {input.reference} \
        {params.name} &&
        echo "`date -R`: Success! Indexes are build." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

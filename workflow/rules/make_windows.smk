
rule make_windows:
    input:
        region="results/regions/{sample}.bed",
        genome="resources/ref_genomes/hg19/genome_hg19.bed",
    output:
        final="results/regions/{sample}_windows_201_100.bed",
    params:
        intervalLen="100",
        windowNum="201",
        rev="",
    log:
        "logs/{sample}_make_windows.log",
    benchmark:
        "logs/{sample}_make_windows.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:  
        """
        (echo "`date -R`: Creating windows..." &&
        workflow/scripts/makewindows.sh \
        {input.region} \
        {input.genome} \
        {output} \
        {params.intervalLen} {params.windowNum} \
        {params.rev} &&
        echo "`date -R`: Success! Windowed file is ready." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """
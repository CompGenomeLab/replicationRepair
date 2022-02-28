
rule make_windows:
    input:
        region="results/regions/{sample}.bed",
        genome="resources/ref_genomes/hg19/genome_hg19.bed",
    output:
        final="results/regions/{sample}_windows_{windowNum}_{interval}.bed",
    params:
        intervalLen="{interval}",
        windowNum="{windowNum}",
        rev="",
    log:
        "logs/rule/analysis/{sample}_{windowNum}_{interval}_make_windows.log",
    benchmark:
        "logs/rule/analysis/{sample}_{windowNum}_{interval}_make_windows.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
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
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """
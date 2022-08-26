
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


rule make_windows_tss_tes:
    input:
        region="resources/ref_genomes/hg19/hg19_ucsc_genes_knownCanonical_stranded.bed",
        genome="resources/ref_genomes/hg19/genome_hg19.bed",
    output:
        tes="results/regions/hg19_ucsc_genes_knownCanonical_tes_windows_2001_1.bed",
        tss="results/regions/hg19_ucsc_genes_knownCanonical_tss_windows_2001_1.bed",
    params:
        intervalLen="1",
        windowNum="2001",
    log:
        "logs/rule/analysis/make_windows_2001_1_tss_tes.log",
    benchmark:
        "logs/rule/analysis/make_windows_2001_1_tss_tes.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Creating tss windows..." &&
        workflow/scripts/makewindows_tes.sh \
        {input.region} \
        {input.genome} \
        {output.tss} \
        {params.intervalLen} {params.windowNum} \
        "" "tss" &&
        echo "`date -R`: Success! Windowed file is ready." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Creating tes windows..." &&
        workflow/scripts/makewindows_tes.sh \
        {input.region} \
        {input.genome} \
        {output.tes} \
        {params.intervalLen} {params.windowNum} \
        "reverse" "tes" &&
        echo "`date -R`: Success! Windowed file is ready." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
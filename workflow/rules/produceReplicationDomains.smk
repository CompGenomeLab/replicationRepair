rule produceReplicationDomains:
    input:  
        norm="results/edu/merge_RT.txt",           
    output:
        norm_loess="results/regions/avg_T_Loess.bed",
        norm="results/regions/repdomains_mean0.5.bed",  
    log:
        "logs/rule/analysis/produceReplicationDomains.log",
    benchmark:
        "logs/rule/analysis/produceReplicationDomains.benchmark.txt",
    conda:
        "../envs/produceReplicationDomains.yaml",
    shell:
        """
        Rscript workflow/scripts/produceRepdomains.R \
        -i {input.norm} \
        -l {output.norm_loess} \
        -o {output.norm} &> {log}
        """

rule produceReplicationDomains_uv:
    input:  
        uv="results/edu/merge_uv_RT.txt",           
    output:
        uv_loess="results/regions/avg_uv_T_Loess.bed", 
        uv="results/regions/repdomains_uv_mean0.5.bed",  
    log:
        "logs/rule/analysis/produceReplicationDomains_uv.log",
    benchmark:
        "logs/rule/analysis/produceReplicationDomains_uv.benchmark.txt",
    conda:
        "../envs/produceReplicationDomains.yaml",
    shell:
        """
        Rscript workflow/scripts/produceRepdomains.R \
        -i {input.uv} \
        -l {output.uv_loess} \
        -o {output.uv} &> {log}
        """
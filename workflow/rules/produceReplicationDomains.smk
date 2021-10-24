rule produceReplicationDomains:
    input:  
        norm="results/edu/merge_RT.txt",  
        uv="results/edu/merge_uv_RT.txt",           
    output:
        norm="results/regions/repdomains_mean0.5.bed",  
        uv="results/regions/repdomains_uv_mean0.5.bed",  
    log:
        "logs/produceReplicationDomains.log",
    benchmark:
        "logs/produceReplicationDomains.benchmark.txt",
    conda:
        "../envs/produceReplicationDomains.yaml",
    shell:
        """
        Rscript workflow/scripts/produceRepdomains.R \
        -i {input.norm} \
        -o {output.norm} &> {log}

        Rscript workflow/scripts/produceRepdomains_uv.R \
        -i {input.uv} \
        -o {output.uv} &>> {log}
        """
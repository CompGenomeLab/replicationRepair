rule bamCompare:
    input:
        early="results/edu/R21061297-EdUrep-UV1-5hrls_combined/R21061297-EdUrep-UV1-5hrls_combined_hg19_sorted_rmdup.bam",
        late="results/edu/R21061297-EdUrep-UV3-5hrls_combined/R21061297-EdUrep-UV3-5hrls_combined_hg19_sorted_rmdup.bam",
    output:
        "results/edu/early_late.bdg",
    log:
        "logs/bamCompare/bamCompare.log",
    benchmark:
        "logs/bamCompare/bamCompare.benchmark.txt",
    conda:
        "../envs/deeptools.yaml"
    shell:  
        """
        (echo "`date -R`: Compare bam files..." &&
        bamCompare \
        -b1 {input.early} \
        -b2 {input.late} \
        --binSize 100000 \
        -o {output} \
        -of "bedgraph" &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

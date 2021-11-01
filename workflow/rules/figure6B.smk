rule figure6B:
    input:  
        noverlap="results/okseq/HeLa_no_overlap.bed",
        overlap="results/okseq/HeLa_intersect2_GM06990_IMR90.bed",       
    output:
        "results/plots/figure6B.pdf",
    log:
        "logs/figure6B.log",
    benchmark:
        "logs/figure6B.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure6B.R \
        --noverlap {input.noverlap} \
        --overlap {input.overlap} \
        --fig6B {output} &> {log}
        """
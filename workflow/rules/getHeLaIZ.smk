rule getHeLaIZ:
    input:  
        hela="results/okseq/HeLa_intersect2_GM06990_IMR90.txt",       
    output:
        noverlap="results/okseq/HeLa_no_overlap.bed",
        overlap="results/okseq/HeLa_intersect2_GM06990_IMR90.bed",
    log:
        "logs/getHeLaIZ.log",
    benchmark:
        "logs/getHeLaIZ.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/getHeLaIZ.R \
        --hela {input.hela} \
        --noverlap {output.noverlap} \
        --overlap {output.overlap} &> {log}
        """
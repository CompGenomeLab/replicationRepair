rule figure_markers:
    input:  
        "results/final/final_reports_markers_iz_repdomains_uv_m0.5_hela_windows_201_100_intergenic.txt",          
    output:
        "results/plots/figure_markers.pdf",
    log:
        "logs/figure_markers.log",
    benchmark:
        "logs/figure_markers.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure_markers.R \
        -i {input} \
        -o {output} &> {log}
        """
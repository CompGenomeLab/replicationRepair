rule figure_markers:
    input:  
        histone="results/final/final_reports_markers_iz_repdomains_uv_m0.5_hela_windows_201_100_intergenic.txt", 
        methyl="results/final/final_reports_methyl_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt",         
    output:
        histone="results/plots/figure_markers.pdf",
        methyl="results/plots/figure_methyl.pdf",
    log:
        "logs/figure_markers.log",
    benchmark:
        "logs/figure_markers.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure_markers.R \
        -i {input.histone} \
        -o {output.histone} &> {log}

        Rscript workflow/scripts/figure_methyl.R \
        -i {input.methyl} \
        -o {output.methyl} &>> {log}
        """
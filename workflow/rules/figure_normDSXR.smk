rule figure_normDSXR:
    input:  
        df="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        df_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",            
    output:
        "results/plots/figure_normDSXR.pdf",
    log:
        "logs/figure_normDSXR.log",
    benchmark:
        "logs/figure_normDSXR.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/normalized_DS_XR.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        -o {output} &> {log}
        """
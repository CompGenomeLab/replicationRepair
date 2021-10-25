rule figure4_5:
    input:  
        df="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        df_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",  
        df_int="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt", 
        df_int_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt",          
    output:
        fig4="results/plots/figure4.pdf",
        fig5="results/plots/figure5.pdf",
    log:
        "logs/figure4_5.log",
    benchmark:
        "logs/figure4_5.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure4.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --intergenic False \
        -o {output.fig4} &> {log}

        Rscript workflow/scripts/figure4.R \
        --df {input.df_int} \
        --df_sim {input.df_int_sim} \
        --intergenic True \
        -o {output.fig5} &>> {log}
        """
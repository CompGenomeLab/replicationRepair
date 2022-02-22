rule figure4C_4D_S10_S11_S12B_S12C:
    input:  
        df="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        df_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",  
        df_int="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt", 
        df_int_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt",          
    output:
        fig4="results/plots/figure4C_4D.pdf",
        fig5="results/plots/figureS11.pdf",
        figS6="results/plots/figureS10.pdf",
        figS7="results/plots/figureS12B_S12C.pdf",
    log:
        "logs/figure4C_4D_S10_S11_S12B_S12C.log",
    benchmark:
        "logs/figure4C_4D_S10_S11_S12B_S12C.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure4C_4D_S10_S11_S12B_S12C.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --phase "late" \
        --intergenic False \
        -o {output.fig4} &> {log}

        Rscript workflow/scripts/figure4C_4D_S10_S11_S12B_S12C.R \
        --df {input.df_int} \
        --df_sim {input.df_int_sim} \
        --phase "late" \
        --intergenic True \
        -o {output.fig5} &>> {log}

        Rscript workflow/scripts/figure4C_4D_S10_S11_S12B_S12C.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --phase "early" \
        --intergenic False \
        -o {output.figS6} &>> {log}

        Rscript workflow/scripts/figure4C_4D_S10_S11_S12B_S12C.R \
        --df {input.df_int} \
        --df_sim {input.df_int_sim} \
        --phase "early" \
        --intergenic True \
        -o {output.figS7} &>> {log}
        """
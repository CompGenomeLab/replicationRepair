rule figure4C_4D_S10_S11_S12B_S12C:
    input:  
        df="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        df_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",  
        df_int="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt", 
        df_int_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt",          
    output:
        fig4=report("results/plots/figure4C_4D.pdf", caption="../report/figure4C_4D.rst", category="Figures"),
        figS11=report("results/plots/figureS11.pdf", caption="../report/figureS11.rst", category="Supplementary Figures"),
        figS10=report("results/plots/figureS10.pdf", caption="../report/figureS10.rst", category="Supplementary Figures"),
        figS12=report("results/plots/figureS12B_S12C.pdf", caption="../report/figureS12.rst", category="Supplementary Figures"),
    log:
        "logs/rule/fig/figure4C_4D_S10_S11_S12B_S12C.log",
    benchmark:
        "logs/rule/fig/figure4C_4D_S10_S11_S12B_S12C.benchmark.txt",
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
        -o {output.figS11} &>> {log}

        Rscript workflow/scripts/figure4C_4D_S10_S11_S12B_S12C.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --phase "early" \
        --intergenic False \
        -o {output.figS10} &>> {log}

        Rscript workflow/scripts/figure4C_4D_S10_S11_S12B_S12C.R \
        --df {input.df_int} \
        --df_sim {input.df_int_sim} \
        --phase "early" \
        --intergenic True \
        -o {output.figS12} &>> {log}
        """
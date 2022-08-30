rule figure4C_4D_S13_S15_S16:
    input:  
        df="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        df_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",  
        df_int="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt", 
        df_int_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt",          
    output:
        fig4=report("results/plots/figure4C_4D.pdf", caption="../report/figure4C_4D.rst", category="Figures"),
        figS13=report("results/plots/figureS13.pdf", caption="../report/figureS13.rst", category="Supplementary Figures"),
        figS15=report("results/plots/figureS15.pdf", caption="../report/figureS15.rst", category="Supplementary Figures"),
        figS16=report("results/plots/figureS16.pdf", caption="../report/figureS16.rst", category="Supplementary Figures"),
        dfs1=report(expand("results/plot_dataframe/figure4_{df}.csv", df=["C1", "C2", "D1", "D2"]), category="Figure Data"),
        dfs2=report(expand("results/plot_dataframe/figureS16_{df}.csv", df=["B1", "B2", "C1", "C2"]), category="Figure Data"),
        dfs3=report(expand("results/plot_dataframe/figure{num}_{df}.csv", num=["S13", "S15"], df=["A1", "A2", "B1", "B2"]),
        category="Figure Data"),
    log:
        "logs/rule/fig/figure4C_4D_S13_S15_S16.log",
    benchmark:
        "logs/rule/fig/figure4C_4D_S13_S15_S16.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure4C_4D_S13_S15_S16.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --phase "late" \
        --intergenic False \
        --data_prefix "results/plot_dataframe/figure4_" \
        -o {output.fig4} &> {log}

        Rscript workflow/scripts/figure4C_4D_S13_S15_S16.R \
        --df {input.df_int} \
        --df_sim {input.df_int_sim} \
        --phase "late" \
        --intergenic True \
        --data_prefix "results/plot_dataframe/figureS15_" \
        -o {output.figS15} &>> {log}

        Rscript workflow/scripts/figure4C_4D_S13_S15_S16.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --phase "early" \
        --intergenic False \
        --data_prefix "results/plot_dataframe/figureS13_" \
        -o {output.figS13} &>> {log}

        Rscript workflow/scripts/figure4C_4D_S13_S15_S16.R \
        --df {input.df_int} \
        --df_sim {input.df_int_sim} \
        --phase "early" \
        --intergenic True \
        --data_prefix "results/plot_dataframe/figureS16_" \
        -o {output.figS16} &>> {log}
        """
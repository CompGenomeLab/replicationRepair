rule figureS11_S12:
    input:  
        df="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        df_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",  
        df_int="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt", 
        df_int_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt",          
    output:
        figS11=report("results/plots/figureS11.pdf", caption="../report/figureS11.rst", category="Supplementary Figures"),
        figS12=report("results/plots/figureS12.pdf", caption="../report/figureS12.rst", category="Supplementary Figures"),
        dfs=report(expand("results/plot_dataframe/figure{num}_{df}.csv", num=["S11", "S12"], df=["A1", "A2", "B1", "B2"]),
        category="Figure Data"),
    log:
        "logs/rule/fig/figureS11_S12.log",
    benchmark:
        "logs/rule/fig/figureS11_S12.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS11_S12.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --intergenic False \
        --data_prefix "results/plot_dataframe/figureS11_" \
        -o {output.figS11} &> {log}

        Rscript workflow/scripts/figureS11_S12.R \
        --df {input.df_int} \
        --df_sim {input.df_int_sim} \
        --intergenic True \
        --data_prefix "results/plot_dataframe/figureS12_" \
        -o {output.figS12} &>> {log}
        """
rule figureS5:
    input:  
        df="results/final/final_reports_hg19_avg_uv_T_Loess.txt", 
        df_sim="results/final/final_reports_sim_hg19_avg_uv_T_Loess.txt",            
    output:
        plot=report("results/plots/figureS5.pdf", caption="../report/figureS5.rst", category="Supplementary Figures"),
        dfs=expand("results/plot_dataframe/figureS5_{df}.csv", df=["A","B"]),
    log:
        "logs/rule/fig/figureS5.log",
    benchmark:
        "logs/rule/fig/figureS5.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS5.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --data_prefix "results/plot_dataframe/figureS5_" \
        -o {output.plot} &> {log}
        """
rule figureS4:
    input:  
        df="results/final/final_reports_hg19_repdomains_uv_mean0.5_windows_201_10000.txt", 
        df_sim="results/final/final_reports_sim_hg19_repdomains_uv_mean0.5_windows_201_10000.txt",            
    output:
        plot=report("results/plots/figureS4.pdf", caption="../report/figureS4.rst", category="Supplementary Figures"),
        dfs="results/plot_dataframe/figureS4.csv"
    log:
        "logs/rule/fig/figureS4.log",
    benchmark:
        "logs/rule/fig/figureS4.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS4.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --data_prefix "results/plot_dataframe/figureS4" \
        -o {output.plot} &> {log}
        """
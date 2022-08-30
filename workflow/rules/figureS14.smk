rule figureS14:
    input:  
        df="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        df_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",            
    output:
        plot=report("results/plots/figureS14.pdf", caption="../report/figureS14.rst", category="Supplementary Figures"),
        dfs=report("results/plot_dataframe/figureS14.csv", category="Figure Data"),
    log:
        "logs/rule/fig/figureS14.log",
    benchmark:
        "logs/rule/fig/figureS14.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS14.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --data_prefix "results/plot_dataframe/figureS14" \
        -o {output.plot} &> {log}
        """
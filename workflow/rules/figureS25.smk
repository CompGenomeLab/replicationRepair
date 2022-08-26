rule figureS25:
    input:  
        df="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        df_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",            
    output:
        plot=report("results/plots/figureS25.pdf", caption="../report/figureS25.rst", category="Supplementary Figures"),
        dfs="results/plot_dataframe/figureS25.csv",
    log:
        "logs/rule/fig/figureS25.log",
    benchmark:
        "logs/rule/fig/figureS25.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS25.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --data_prefix "results/plot_dataframe/figureS25" \
        -o {output.plot} &> {log}
        """
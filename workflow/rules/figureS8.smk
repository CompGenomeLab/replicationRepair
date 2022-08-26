rule figureS8:
    input:  
        df="results/plot_dataframe/figure3_B.csv",            
    output:
        plot=report("results/plots/figureS8.pdf", caption="../report/figureS8.rst", category="Supplementary Figures"),
        dfs="results/plot_dataframe/figureS8.csv",
    log:
        "logs/rule/fig/figureS8.log",
    benchmark:
        "logs/rule/fig/figureS8.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS8.R \
        --df {input.df} \
        --data_prefix "results/plot_dataframe/figureS8" \
        -o {output.plot} &> {log}
        """
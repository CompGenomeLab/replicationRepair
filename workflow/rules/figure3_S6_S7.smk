rule figure3_S6_S7:
    input:  
        df="results/final/final_reports_hg19_chromhmm_hela_repdomains_uv_mean0.5.txt", 
        df_sim="results/final/final_reports_sim_hg19_chromhmm_hela_repdomains_uv_mean0.5.txt",            
    output:
        p3=report("results/plots/figure3.pdf", caption="../report/figure3.rst", category="Figures"),
        S6=report("results/plots/figureS6.pdf", caption="../report/figureS6.rst", category="Supplementary Figures"),
        S7=report("results/plots/figureS7.pdf", caption="../report/figureS7.rst", category="Supplementary Figures"),
        dfs=expand("results/plot_dataframe/figure{num}_{df}.csv", num=["3", "S6", "S7"], df=["A", "B"]),
    log:
        "logs/rule/fig/figure3_S6_S7.log",
    benchmark:
        "logs/rule/fig/figure3_S6_S7.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure3_S6.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --dtype "CPD" \
        --data_prefix "results/plot_dataframe/figure3_" \
        -o {output.p3} &> {log}

        Rscript workflow/scripts/figure3_S6.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --dtype "64_PP" \
        --data_prefix "results/plot_dataframe/figureS6_" \
        -o {output.S6} &>> {log}

        Rscript workflow/scripts/figureS7.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --data_prefix "results/plot_dataframe/figureS7_" \
        -o {output.S7} &>> {log}
        """
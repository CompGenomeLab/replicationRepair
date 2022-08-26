rule figureS9_S10:
    input:  
        real="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",   
    output:
        figS9=report("results/plots/figureS9.pdf", caption="../report/figureS9.rst", category="Supplementary Figures"),
        figS10=report("results/plots/figureS10.pdf", caption="../report/figureS10.rst", category="Supplementary Figures"),
        dfs=expand("results/plot_dataframe/figure{num}_{df}.csv", num=["S9", "S10"], df=["A1", "A2", "B1", "B2"]),
    log:
        "logs/rule/fig/figureS9_S10.log",
    benchmark:
        "logs/rule/fig/figureS9_S10.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS9_S10.R \
        --real {input.real} \
        --sim {input.sim} \
        --prod CPD \
        --data_prefix "results/plot_dataframe/figureS9_" \
        -o {output.figS9} &> {log}

        Rscript workflow/scripts/figureS9_S10.R \
        --real {input.real} \
        --sim {input.sim} \
        --prod 64_PP \
        --data_prefix "results/plot_dataframe/figureS10_" \
        -o {output.figS10} &>> {log}
        """
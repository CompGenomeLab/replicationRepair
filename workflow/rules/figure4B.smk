
rule figure4B:
    input:  
        real="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",  
    output:
        plot=report("results/plots/figure4B.pdf", caption="../report/figure4B.rst", category="Figures"), 
        dfs=expand("results/plot_dataframe/figure4_{df}.csv", df=["B1", "B2", "B3", "B4"]),
    log:
        "logs/rule/fig/figure4B.log",
    benchmark:
        "logs/rule/fig/figure4B.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure4B.R \
        --real {input.real} \
        --sim {input.sim} \
        --prod CPD \
        --data_prefix "results/plot_dataframe/figure4_" \
        -o {output.plot} &> {log}
        """
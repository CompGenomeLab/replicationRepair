rule figure2:
    input:  
        windowed="results/final/final_reports_hg19_repdomains_uv_mean0.5_windows_201_10000.txt", 
        windowed_sim="results/final/final_reports_sim_hg19_repdomains_uv_mean0.5_windows_201_10000.txt",        
        noW="results/final/final_reports_hg19_repdomains_uv_mean0.5.txt",            
        noW_sim="results/final/final_reports_sim_hg19_repdomains_uv_mean0.5.txt",                
    output:
        plot=report("results/plots/figure2.pdf", caption="../report/figure2.rst", category="Figures"),
        dfs=report(expand("results/plot_dataframe/figure2_{df}.csv", df=["A1", "A2", "B"]),
        category="Figure Data"),
    log:
        "logs/rule/fig/figure2.log",
    benchmark:
        "logs/rule/fig/figure2.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure2.R \
        --windowed {input.windowed} \
        --windowed_sim {input.windowed_sim} \
        --noW {input.noW} \
        --noW_sim {input.noW_sim} \
        --data_prefix "results/plot_dataframe/figure2_" \
        -o {output.plot} &> {log}
        """


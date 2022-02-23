rule figureS4:
    input:  
        df="results/final/final_reports_hg19_repdomains_uv_mean0.5_windows_201_10000.txt", 
        df_sim="results/final/final_reports_sim_hg19_repdomains_uv_mean0.5_windows_201_10000.txt",            
    output:
        "results/plots/figureS4.pdf",
    log:
        "logs/figureS4.log",
    benchmark:
        "logs/figureS4.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS4.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        -o {output} &> {log}
        """
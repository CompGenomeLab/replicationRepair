rule figureS18:
    input:  
        df="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        df_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",            
    output:
        report("results/plots/figureS18.pdf", caption="../report/figureS18.rst", category="Supplementary Figures"),
    log:
        "logs/figureS18.log",
    benchmark:
        "logs/figureS18.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS18.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        -o {output} &> {log}
        """
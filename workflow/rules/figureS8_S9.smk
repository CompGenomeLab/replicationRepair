rule figureS8_S9:
    input:  
        df="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        df_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",  
        df_int="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt", 
        df_int_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt",          
    output:
        figS8=report("results/plots/figureS8.pdf", caption="../report/figureS8.rst", category="Supplementary Figures"),
        figS9=report("results/plots/figureS9.pdf", caption="../report/figureS9.rst", category="Supplementary Figures"),
    log:
        "logs/rule/fig/figureS8_S9.log",
    benchmark:
        "logs/rule/fig/figureS8_S9.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS8_S9.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --intergenic False \
        -o {output.figS8} &> {log}

        Rscript workflow/scripts/figureS8_S9.R \
        --df {input.df_int} \
        --df_sim {input.df_int_sim} \
        --intergenic True \
        -o {output.figS9} &>> {log}
        """
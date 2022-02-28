rule figure3_S5:
    input:  
        df="results/final/final_reports_hg19_chromhmm_hela_repdomains_uv_mean0.5.txt", 
        df_sim="results/final/final_reports_sim_hg19_chromhmm_hela_repdomains_uv_mean0.5.txt",            
    output:
        CPD=report("results/plots/figure3.pdf", caption="../report/figure3.rst", category="Figures"),
        PP64=report("results/plots/figureS5.pdf", caption="../report/figureS5.rst", category="Supplementary Figures"),
    log:
        "logs/rule/fig/figure3.log",
    benchmark:
        "logs/rule/fig/figure3.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure3_S5.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --dtype "CPD" \
        -o {output.CPD} &> {log}

        Rscript workflow/scripts/figure3_S5.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --dtype "64_PP" \
        -o {output.PP64} &>> {log}
        """
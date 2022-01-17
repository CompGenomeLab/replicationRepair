rule figure3:
    input:  
        df="results/final/final_reports_hg19_chromhmm_hela_repdomains_uv_mean0.5.txt", 
        df_sim="results/final/final_reports_sim_hg19_chromhmm_hela_repdomains_uv_mean0.5.txt",            
    output:
        CPD="results/plots/figure3.pdf",
        PP64="results/plots/figureS4.pdf",
        CPD_ttest="results/plots/figure3_ttest.pdf",
        PP64_ttest="results/plots/figureS4_ttest.pdf",
    log:
        "logs/figure3.log",
    benchmark:
        "logs/figure3.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure3.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --dtype "CPD" \
        --o1 {output.CPD} \
        --o2 {output.CPD_ttest} &> {log}

        Rscript workflow/scripts/figure3.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --dtype "64_PP" \
        --o1 {output.PP64} \
        --o2 {output.PP64_ttest} &>> {log}
        """
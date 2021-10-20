rule figure3:
    input:  
        df="results/final/final_reports_hg19_chromhmm_repdomains_hela_windows_chr.txt", 
        df_sim="results/final/final_reports_sim_hg19_chromhmm_repdomains_hela_windows_chr.txt",            
    output:
        CPD="results/plots/figure3.pdf",
        PP64="results/plots/supplementary_figure_ChromHMM_64PP.pdf",
    log:
        "logs/figure3.log",
    benchmark:
        "logs/figure3.benchmark.txt",
    conda:
        "../envs/figure1_2_3.yaml",
    shell:
        """
        Rscript workflow/scripts/figure3.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --dtype "CPD" \
        -o {output.CPD} &> {log}

        Rscript workflow/scripts/figure3.R \
        --df {input.df} \
        --df_sim {input.df_sim} \
        --dtype "64_PP" \
        -o {output.PP64} &>> {log}
        """
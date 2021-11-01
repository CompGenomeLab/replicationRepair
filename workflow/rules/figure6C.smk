rule figure6C:
    input:  
        "results/intergenic_mutation/melanoma/melanoma_target_mut_iz_hela_to_gm_imr_repdomains_uv_mean0.5_with_scores_windows_201_100_combined_rpkm.txt",   
    output:
        "results/plots/figure6C.pdf",
    log:
        "logs/figure6C.log",
    benchmark:
        "logs/figure6C.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure6C.R \
        --hela {input} \
        --fig6C {output} &> {log}
        """
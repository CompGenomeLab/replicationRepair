rule figure6D:
    input:  
        overlap="results/intergenic_mutation/melanoma/melanoma_target_mut_iz_hela_to_gm_imr_repdomains_uv_mean0.5_windows_201_100_combined_rpkm.txt",
        noverlap="results/intergenic_mutation/melanoma/melanoma_target_mut_iz_hela_no_overlap_repdomains_uv_mean0.5_windows_201_100_combined_rpkm.txt",  
        overlap_TC="results/regions/iz_hela_to_gm_imr_repdomains_uv_mean0.5_windows_201_100_counts.txt",
        noverlap_TC="results/regions/iz_hela_no_overlap_repdomains_uv_mean0.5_windows_201_100_counts.txt",        
    output:
        "results/plots/figure6D.pdf",
    log:
        "logs/figure6D.log",
    benchmark:
        "logs/figure6D.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure6D.R \
        --noverlap {input.noverlap} \
        --overlap {input.overlap} \
        --noverlap_TC {input.noverlap_TC} \
        --overlap_TC {input.overlap_TC} \
        --fig6D {output} &> {log}
        """
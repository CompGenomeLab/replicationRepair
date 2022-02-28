rule figure5B_5C_5D:
    input:  
        noverlap="results/okseq/HeLa_no_overlap.bed",
        overlap="results/okseq/HeLa_intersect2_GM06990_IMR90.bed",  
        hela="results/intergenic_mutation/melanoma/melanoma_target_mut_iz_hela_to_gm_imr_repdomains_uv_mean0.5_with_scores_windows_201_100_combined_rpkm.txt",
        overlap_mut="results/intergenic_mutation/melanoma/melanoma_target_mut_iz_hela_to_gm_imr_repdomains_uv_mean0.5_windows_201_100_combined_rpkm.txt",
        noverlap_mut="results/intergenic_mutation/melanoma/melanoma_target_mut_iz_hela_no_overlap_repdomains_uv_mean0.5_windows_201_100_combined_rpkm.txt",  
        overlap_TC="results/regions/iz_hela_to_gm_imr_repdomains_uv_mean0.5_windows_201_100_counts.txt",
        noverlap_TC="results/regions/iz_hela_no_overlap_repdomains_uv_mean0.5_windows_201_100_counts.txt",  
    output:
        report("results/plots/figure5B_5C_5D.pdf", caption="../report/figure5B_5C_5D.rst", category="Figures"),
    log:
        "logs/rule/fig/figure5B_5C_5D.log",
    benchmark:
        "logs/rule/fig/figure5B_5C_5D.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure5B_5C_5D.R \
        --noverlap {input.noverlap} \
        --overlap {input.overlap} \
        --hela {input.hela} \
        --noverlap_mut {input.noverlap_mut} \
        --overlap_mut {input.overlap_mut} \
        --noverlap_TC {input.noverlap_TC} \
        --overlap_TC {input.overlap_TC} \
        --fig5 {output} &> {log}
        """
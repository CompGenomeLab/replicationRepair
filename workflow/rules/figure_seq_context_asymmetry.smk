rule figure_seq_context_asymmetry:
    input:  
        "results/regions/iz_hela_repdomains_uv_mean0.5_windows_201_100_kmer_counts.txt", 
    output:
        comb="results/plots/figure4_seq_asymmetry.pdf",
        table="results/table/seq_asymmetry_exp_obs.csv",
    log:
        "logs/figure_seq_context_asymmetry.log",
    benchmark:
        "logs/figure_seq_context_asymmetry.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/nucleotide_asymmetry.R \
        -i {input} \
        --comb {output.comb} \
        --exp_obs {output.table} &> {log}
        """
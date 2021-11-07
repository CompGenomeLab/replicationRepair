rule figure_seq_context_asymmetry:
    input:  
        "results/regions/iz_hela_repdomains_uv_mean0.5_windows_201_100_kmer_counts.txt", 
    output:
        o1="results/plots/mer1.pdf",
        o2="results/plots/mer2.pdf",
        o3="results/plots/mer3.pdf",
        o4="results/plots/mer4.pdf",
        o5="results/plots/mer5.pdf",
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
        --o1 {output.o1} \
        --o2 {output.o2} \
        --o3 {output.o3} \
        --o4 {output.o4} \
        --o5 {output.o5} &> {log}
        """
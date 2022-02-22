rule figure4E_S13_S14_S15_S16_S17:
    input:  
        "results/regions/iz_hela_repdomains_uv_mean0.5_windows_201_100_kmer_counts.txt", 
    output:
        comb="results/plots/figure4E.pdf",
        table="results/table/seq_asymmetry_exp_obs.csv",
        o1="results/plots/figureS13.pdf",
        o2="results/plots/figureS14.pdf",
        o3="results/plots/figureS15.pdf",
        o4="results/plots/figureS16.pdf",
        o5="results/plots/figureS17.pdf",
    log:
        "logs/figure4E_S13_S14_S15_S16_S17.log",
    benchmark:
        "logs/figure4E_S13_S14_S15_S16_S17.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure4E_S13_S14_S15_S16_S17.R \
        -i {input} \
        --comb {output.comb} \
        --o1 {output.o1} \
        --o2 {output.o2} \
        --o3 {output.o3} \
        --o4 {output.o4} \
        --o5 {output.o5} \
        --exp_obs {output.table} &> {log}
        """
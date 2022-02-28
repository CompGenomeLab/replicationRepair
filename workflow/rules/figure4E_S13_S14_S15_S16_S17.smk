rule figure4E_S13_S14_S15_S16_S17:
    input:  
        "results/regions/iz_hela_repdomains_uv_mean0.5_windows_201_100_kmer_counts.txt", 
    output:
        comb=report("results/plots/figure4E.pdf", caption="../report/figure4E.rst", category="Figures"),
        table=report("results/table/seq_asymmetry_exp_obs.csv", caption="../report/table.rst", category="Table"),
        o1=report("results/plots/figureS13.pdf", caption="../report/figureS13.rst", category="Supplementary Figures"),
        o2=report("results/plots/figureS14.pdf", caption="../report/figureS14.rst", category="Supplementary Figures"),
        o3=report("results/plots/figureS15.pdf", caption="../report/figureS15.rst", category="Supplementary Figures"),
        o4=report("results/plots/figureS16.pdf", caption="../report/figureS16.rst", category="Supplementary Figures"),
        o5=report("results/plots/figureS17.pdf", caption="../report/figureS17.rst", category="Supplementary Figures"),
    log:
        "logs/rule/fig/figure4E_S13_S14_S15_S16_S17.log",
    benchmark:
        "logs/rule/fig/figure4E_S13_S14_S15_S16_S17.benchmark.txt",
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
rule figure4E_4F_S17toS24:
    input:  
        "results/regions/iz_hela_repdomains_uv_mean0.5_windows_201_100_kmer_counts.txt", 
    output:
        oT=report("results/plots/figure4E.pdf", caption="../report/figure4E.rst", category="Figures"),
        oTT=report("results/plots/figure4F.pdf", caption="../report/figure4F.rst", category="Figures"),
        oC=report("results/plots/figureS23.pdf", caption="../report/figureS23.rst", category="Supplementary Figures"),
        oCC=report("results/plots/figureS24.pdf", caption="../report/figureS24.rst", category="Supplementary Figures"),
        oTC=report("results/plots/figureS22.pdf", caption="../report/figureS22.rst", category="Supplementary Figures"),
        o1=report("results/plots/figureS17.pdf", caption="../report/figureS17.rst", category="Supplementary Figures"),
        o2=report("results/plots/figureS18.pdf", caption="../report/figureS18.rst", category="Supplementary Figures"),
        o3=report("results/plots/figureS19.pdf", caption="../report/figureS19.rst", category="Supplementary Figures"),
        o4=report("results/plots/figureS20.pdf", caption="../report/figureS20.rst", category="Supplementary Figures"),
        o5=report("results/plots/figureS21.pdf", caption="../report/figureS21.rst", category="Supplementary Figures"),
        dfs=report(expand("results/plot_dataframe/figure{num}.csv", num=["S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "4_E", "4_F"]),
        category="Figure Data"),
    log:
        "logs/rule/fig/figure4E_4F_S17toS24.log",
    benchmark:
        "logs/rule/fig/figure4E_4F_S17toS24.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure4E_4F_S17toS24.R \
        -i {input} \
        --data_prefix "results/plot_dataframe/figure" \
        --oT {output.oT} \
        --oTT {output.oTT} \
        --oC {output.oC} \
        --oCC {output.oCC} \
        --oTC {output.oTC} \
        --o1 {output.o1} \
        --o2 {output.o2} \
        --o3 {output.o3} \
        --o4 {output.o4} \
        --o5 {output.o5} &> {log}
        """
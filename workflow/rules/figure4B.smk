
rule figure4B:
    input:  
        real="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",  
    output:
        report("results/plots/figure4B.pdf", caption="../report/figure4B.rst", category="Figures"), 
    log:
        "logs/figure4B.log",
    benchmark:
        "logs/figure4B.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure4B.R \
        --real {input.real} \
        --sim {input.sim} \
        --prod CPD \
        -o {output} &> {log}
        """
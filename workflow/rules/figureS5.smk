rule figureS5:
    input:  
        real="results/final/final_reports_hg19_iz_hela_windows_201_100.txt", 
        sim="results/final/final_reports_sim_hg19_iz_hela_windows_201_100.txt",          
    output:
        "results/plots/figureS5.pdf",
    log:
        "logs/figureS5.log",
    benchmark:
        "logs/figureS5.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS5.R \
        --real {input.real} \
        --sim {input.sim} \
        -o {output} &> {log}
        """
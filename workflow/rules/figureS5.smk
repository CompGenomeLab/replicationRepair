rule figureS5:
    input:  
        real="results/final/final_reports_hg19_iz_hela_windows_201_100.txt", 
        sim="results/final/final_reports_sim_hg19_iz_hela_windows_201_100.txt",  
        int_real="results/final/final_reports_hg19_iz_hela_windows_201_100_intergenic.txt", 
        int_sim="results/final/final_reports_sim_hg19_iz_hela_windows_201_100_intergenic.txt",        
    output:
        figS5="results/plots/figureS5.pdf",
        int_figS5="results/plots/figureS5_intergenic.pdf",
        figS5_64="results/plots/figureS5_64.pdf",
        int_figS5_64="results/plots/figureS5_intergenic_64.pdf",        
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
        --prod CPD \
        -o {output.figS5} &> {log}

        Rscript workflow/scripts/figureS5.R \
        --real {input.int_real} \
        --sim {input.int_sim} \
        --prod CPD \
        -o {output.int_figS5} &>> {log}

        Rscript workflow/scripts/figureS5.R \
        --real {input.real} \
        --sim {input.sim} \
        --prod 64_PP \
        -o {output.figS5_64} &>> {log}

        Rscript workflow/scripts/figureS5.R \
        --real {input.int_real} \
        --sim {input.int_sim} \
        --prod 64_PP \
        -o {output.int_figS5_64} &>> {log}
        """
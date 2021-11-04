rule figureS5_repdomains:
    input:  
        real="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",  
        int_real="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt", 
        int_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt",   
    output:
        figS5="results/plots/figureS5_repdomains.pdf",
        int_figS5="results/plots/figureS5_repdomains_intergenic.pdf",
        figS5_64="results/plots/figureS5_repdomains_64.pdf",
        int_figS5_64="results/plots/figureS5_repdomains_intergenic_64.pdf",        
    log:
        "logs/figureS5_repdomains.log",
    benchmark:
        "logs/figureS5_repdomains.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS5_repdomains.R \
        --real {input.real} \
        --sim {input.sim} \
        --prod CPD \
        -o {output.figS5} &> {log}

        Rscript workflow/scripts/figureS5_repdomains.R \
        --real {input.int_real} \
        --sim {input.int_sim} \
        --prod CPD \
        -o {output.int_figS5} &>> {log}

        Rscript workflow/scripts/figureS5_repdomains.R \
        --real {input.real} \
        --sim {input.sim} \
        --prod 64_PP \
        -o {output.figS5_64} &>> {log}

        Rscript workflow/scripts/figureS5_repdomains.R \
        --real {input.int_real} \
        --sim {input.int_sim} \
        --prod 64_PP \
        -o {output.int_figS5_64} &>> {log}
        """
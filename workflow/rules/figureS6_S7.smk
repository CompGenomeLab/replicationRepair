rule figureS6_S7:
    input:  
        real="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",  
        int_real="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt", 
        int_sim="results/final/final_reports_sim_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100_intergenic.txt",   
    output:
        figS6=report("results/plots/figureS6.pdf", caption="../report/figureS6.rst", category="Supplementary Figures"),
        int_figS6="results/plots/figureS6_intergenic.pdf",
        figS7=report("results/plots/figureS7.pdf", caption="../report/figureS7.rst", category="Supplementary Figures"),
        int_figS7="results/plots/figureS7_intergenic.pdf",        
    log:
        "logs/figureS6_S7.log",
    benchmark:
        "logs/figureS6_S7.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS6_S7.R \
        --real {input.real} \
        --sim {input.sim} \
        --prod CPD \
        -o {output.figS6} &> {log}

        Rscript workflow/scripts/figureS6_S7.R \
        --real {input.int_real} \
        --sim {input.int_sim} \
        --prod CPD \
        -o {output.int_figS6} &>> {log}

        Rscript workflow/scripts/figureS6_S7.R \
        --real {input.real} \
        --sim {input.sim} \
        --prod 64_PP \
        -o {output.figS7} &>> {log}

        Rscript workflow/scripts/figureS6_S7.R \
        --real {input.int_real} \
        --sim {input.int_sim} \
        --prod 64_PP \
        -o {output.int_figS7} &>> {log}
        """
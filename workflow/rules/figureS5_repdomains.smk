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

rule figureS5_repdomains_kmer:
    input:  
        real="results/final/final_reports_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        mer2="results/final/final_reports_sim_2_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        mer3="results/final/final_reports_sim_3_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",  
        mer4="results/final/final_reports_sim_4_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",
        mer5="results/final/final_reports_sim_5_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt", 
        IZ="results/final/final_reports_sim_IZ_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",
        rmTTTT="results/final/final_reports_sim_rmTTTT_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",
        IZ_rmTTTT="results/final/final_reports_sim_IZ_rmTTTT_hg19_iz_hela_repdomains_uv_mean0.5_windows_201_100.txt",
    output:
        mer2="results/plots/figureS5_repdomains_2.pdf",
        mer3="results/plots/figureS5_repdomains_3.pdf",
        mer4="results/plots/figureS5_repdomains_4.pdf",
        mer5="results/plots/figureS5_repdomains_5.pdf",     
        IZ="results/plots/figureS5_repdomains_IZ.pdf",
        rmTTTT="results/plots/figureS5_repdomains_rmTTTT.pdf",    
        IZ_rmTTTT="results/plots/figureS5_repdomains_IZ_rmTTTT.pdf",
    log:
        "logs/figureS5_repdomains_kmer.log",
    benchmark:
        "logs/figureS5_repdomains_kmer.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS5_repdomains.R \
        --real {input.real} \
        --sim {input.mer2} \
        --prod CPD \
        -o {output.mer2} &> {log}

        Rscript workflow/scripts/figureS5_repdomains.R \
        --real {input.real} \
        --sim {input.mer3} \
        --prod CPD \
        -o {output.mer3} &>> {log}

        Rscript workflow/scripts/figureS5_repdomains.R \
        --real {input.real} \
        --sim {input.mer4} \
        --prod CPD \
        -o {output.mer4} &>> {log}

        Rscript workflow/scripts/figureS5_repdomains.R \
        --real {input.real} \
        --sim {input.mer5} \
        --prod CPD \
        -o {output.mer5} &>> {log}

        Rscript workflow/scripts/figureS5_repdomains.R \
        --real {input.real} \
        --sim {input.IZ} \
        --prod CPD \
        -o {output.IZ} &>> {log}

        Rscript workflow/scripts/figureS5_repdomains.R \
        --real {input.real} \
        --sim {input.rmTTTT} \
        --prod CPD \
        -o {output.rmTTTT} &>> {log}
        
        Rscript workflow/scripts/figureS5_repdomains.R \
        --real {input.real} \
        --sim {input.IZ_rmTTTT} \
        --prod CPD \
        -o {output.IZ_rmTTTT} &>> {log}
        """
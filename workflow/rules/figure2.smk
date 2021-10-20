rule figure2:
    input:  
        windowed="results/final/final_reports_hg19_repdomains_uv_mean0.5_windows_201_10000.txt", 
        windowed_sim="results/final/final_reports_sim_hg19_repdomains_uv_mean0.5_windows_201_10000.txt",        
        noW="results/final/final_reports_noWindows_hg19_repdomains_hela.txt",            
        noW_sim="results/final/final_reports_noWindows_sim_hg19_repdomains_hela.txt",                
    output:
        "results/plots/figure2.pdf",
    log:
        "logs/figure2.log",
    benchmark:
        "logs/figure2.benchmark.txt",
    conda:
        "../envs/figure1_2_3.yaml",
    shell:
        """
        Rscript workflow/scripts/figure2.R \
        --windowed {input.windowed} \
        --windowed_sim {input.windowed_sim} \
        --noW {input.noW} \
        --noW_sim {input.noW_sim} \
        -o {output} &> {log}
        """


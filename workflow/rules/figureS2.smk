rule figureS2:
    input:  
        dDSCPD="results/final/XT4H1-CT5H2-HelaD3-5R2h2_combined_cutadapt_sorted_10_dinucleotideTable.txt",           
        dXRCPD="results/final/R19026421-2019XR1-Hela35X3_combined_R1_cutadapt_sorted_26_dinucleotideTable.txt",      
        lenCPD="results/final/R19026421-2019XR1-Hela35X3_combined_R1_cutadapt_length_distribution.txt",            
        tss="results/final/final_reports_hg19_tss.txt",            
        tes="results/final/final_reports_hg19_tes.txt",                
    output:
        "results/plots/figureS2.pdf",
    log:
        "logs/figureS2.log",
    benchmark:
        "logs/figureS2.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS2.R \
        --dinuc_ds_cpd_120 {input.dDSCPD} \
        --dinuc_xr_cpd_120 {input.dXRCPD} \
        --len_xr_cpd_120 {input.lenCPD} \
        --tss {input.tss} \
        --tes {input.tes} \
        -o {output} \
        --log {log} &> {log}
        """


rule figure1:
    input:  
        dDS64="resources/samples/dinuc/HDL64A5_ACAGTG_cutadapt_sorted_10_dinucleotideTable.txt", 
        dDSCPD="resources/samples/dinuc/HDLCA12_CTTGTA_cutadapt_sorted_10_dinucleotideTable.txt",        
        dXR64="resources/samples/dinuc/HXL64A3_TTAGGC_cutadapt_sorted_26_dinucleotideTable.txt",            
        dXRCPD="resources/samples/dinuc/HXLCA6_GCCAAT_cutadapt_sorted_26_dinucleotideTable.txt",            
        len64="resources/samples/dinuc/HXL64A3_TTAGGC_cutadapt_length_distribution.txt",            
        lenCPD="resources/samples/dinuc/HXLCA6_GCCAAT_cutadapt_length_distribution.txt",            
        tss="results/final/final_reports_hg19_tss.txt",            
        tes="results/final/final_reports_hg19_tes.txt",                
    output:
        plot=report("results/plots/figure1.pdf", caption="../report/figure1.rst", category="Figures"),
        dfs=expand("results/plot_dataframe/figure1_{df}.csv", df=["B1", "B2", "C1", "C2", "C3", "C4", "D"]),
    log:
        "logs/rule/fig/figure1.log",
    benchmark:
        "logs/rule/fig/figure1.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure1.R \
        --dinuc_ds_64_12 {input.dDS64} \
        --dinuc_ds_cpd_12 {input.dDSCPD} \
        --dinuc_xr_64_12 {input.dXR64} \
        --dinuc_xr_cpd_12 {input.dXRCPD} \
        --len_xr_64_12 {input.len64} \
        --len_xr_cpd_12 {input.lenCPD} \
        --tss {input.tss} \
        --tes {input.tes} \
        --data_prefix "results/plot_dataframe/figure1_" \
        -o {output.plot} \
        --log {log} 
        """


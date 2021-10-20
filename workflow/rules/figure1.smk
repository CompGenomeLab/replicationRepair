rule figure1:
    input:  
        dDS64="results/final/HDL64A5_ACAGTG_cutadapt_sorted_10_dinucleotideTable.txt", 
        dDSCPD="results/final/HDLCA12_CTTGTA_cutadapt_sorted_10_dinucleotideTable.txt",        
        dXR64="results/final/HXL64A3_TTAGGC_cutadapt_sorted_26_dinucleotideTable.txt",            
        dXRCPD="results/final/HXLCA6_GCCAAT_cutadapt_sorted_26_dinucleotideTable.txt",            
        len64="results/final/HXL64A3_TTAGGC_cutadapt_length_distribution.txt",            
        lenCPD="results/final/HXLCA6_GCCAAT_cutadapt_length_distribution.txt",            
        tss="results/final/final_reports_hg19_tss.txt",            
        tes="results/final/final_reports_hg19_tes.txt",                
    output:
        B="results/plots/figure1B.pdf",
        C="results/plots/figure1C.pdf",
        D="results/plots/figure1D.pdf",
    log:
        "logs/figure1.log",
    benchmark:
        "logs/figure1.benchmark.txt",
    conda:
        "../envs/figure1_2_3.yaml",
    shell:
        """
        Rscript workflow/scripts/figure1_TSS_TES.R \
        --dinuc_ds_64_12 {input.dDS64} \
        --dinuc_ds_cpd_12 {input.dDSCPD} \
        --dinuc_xr_64_12 {input.dXR64} \
        --dinuc_xr_cpd_12 {input.dXRCPD} \
        --len_xr_64_12 {input.len64} \
        --len_xr_cpd_12 {input.lenCPD} \
        --tss {input.tss} \
        --tes {input.tes} \
        -o1 {output.B} \
        -o2 {output.C} \
        -o3 {output.D} \
        --log {log} &> {log}
        """


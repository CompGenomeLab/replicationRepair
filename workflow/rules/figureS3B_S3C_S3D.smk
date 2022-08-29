rule figureS3B_S3C_S3D:
    input:  
        dDSCPD="resources/samples/dinuc/XT4H1-CT5H2-HelaD3-5R2h2_combined_hg19_sorted_dinucleotideTable.txt",           
        dXRCPD="resources/samples/dinuc/R19026421-2019XR1-Hela35X3_combined_R1_hg19_sorted_dinucleotideTable.txt",      
        lenCPD="resources/samples/dinuc/R19026421-2019XR1-Hela35X3_combined_R1_hg19_length_distribution.txt",            
        tss="results/final/final_reports_hg19_tss.txt",            
        tes="results/final/final_reports_hg19_tes.txt",                
    output:
        plot=report("results/plots/figureS3B_S3C_S3D.pdf", caption="../report/figureS3.rst", category="Supplementary Figures"),
        dfs=expand("results/plot_dataframe/figureS3_{df}.csv", df=["B", "C1", "C2", "D"]),
    log:
        "logs/rule/fig/figureS3B_S3C_S3D.log",
    benchmark:
        "logs/rule/fig/figureS3B_S3C_S3D.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figureS3B_S3C_S3D.R \
        --dinuc_ds_cpd_120 {input.dDSCPD} \
        --dinuc_xr_cpd_120 {input.dXRCPD} \
        --len_xr_cpd_120 {input.lenCPD} \
        --tss {input.tss} \
        --tes {input.tes} \
        --data_prefix "results/plot_dataframe/figureS3_" \
        -o {output.plot} \
        --log {log}
        """


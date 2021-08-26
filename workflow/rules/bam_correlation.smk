rule bam_correlation_pe_edu:
    input:
        "results/edu/R21061297-EdUrep-2hrls_combined/R21061297-EdUrep-2hrls_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21061297-EdUrep-4hrls_combined/R21061297-EdUrep-4hrls_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21061297-EdUrep-UV1-5hrls_combined/R21061297-EdUrep-UV1-5hrls_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21061297-EdUrep-UV3-5hrls_combined/R21061297-EdUrep-UV3-5hrls_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21071354-EdUrep2-2hrls2_combined/R21071354-EdUrep2-2hrls2_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21071354-EdUrep2-4hrls2_combined/R21071354-EdUrep2-4hrls2_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21071354-EdUrep2-UV1-5hrls2_combined/R21071354-EdUrep2-UV1-5hrls2_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21071354-EdUrep2-UV3-5hrls2_combined/R21071354-EdUrep2-UV3-5hrls2_combined_hg19_sorted_rmdup.bam",
    output:
        out="results/edu/readCounts.npz",
        raw_out="results/edu/readCounts.tab",
    log:
        "logs/bam_correlation_pe_edu.log",
    benchmark:
        "logs/bam_correlation_pe_edu.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: MultiBam summary..." &&
        multiBamSummary bins \
        --bamfiles {input} \
        --minMappingQuality 20 \
        --labels early_rep1 late_rep1 early_uv_rep1 late_uv_rep1 early_rep2 late_rep2 early_uv_rep2 late_uv_rep2 \
        -out {output.out} --outRawCounts {output.raw_out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule bam_correlation_pe_edu_early:
    input:
        "results/edu/R21061297-EdUrep-2hrls_combined/R21061297-EdUrep-2hrls_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21061297-EdUrep-UV1-5hrls_combined/R21061297-EdUrep-UV1-5hrls_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21071354-EdUrep2-2hrls2_combined/R21071354-EdUrep2-2hrls2_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21071354-EdUrep2-UV1-5hrls2_combined/R21071354-EdUrep2-UV1-5hrls2_combined_hg19_sorted_rmdup.bam",
    output:
        out="results/edu/readCounts_early.npz",
        raw_out="results/edu/readCounts_early.tab",
    log:
        "logs/bam_correlation_pe_edu_early.log",
    benchmark:
        "logs/bam_correlation_pe_edu_early.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: MultiBam summary..." &&
        multiBamSummary bins \
        --bamfiles {input} \
        --minMappingQuality 20 \
        --labels early_rep1 early_uv_rep1 early_rep2 early_uv_rep2 \
        -out {output.out} --outRawCounts {output.raw_out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule bam_correlation_pe_edu_late:
    input:
        "results/edu/R21061297-EdUrep-4hrls_combined/R21061297-EdUrep-4hrls_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21061297-EdUrep-UV3-5hrls_combined/R21061297-EdUrep-UV3-5hrls_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21071354-EdUrep2-4hrls2_combined/R21071354-EdUrep2-4hrls2_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21071354-EdUrep2-UV3-5hrls2_combined/R21071354-EdUrep2-UV3-5hrls2_combined_hg19_sorted_rmdup.bam",
    output:
        out="results/edu/readCounts_late.npz",
        raw_out="results/edu/readCounts_late.tab",
    log:
        "logs/bam_correlation_pe_edu_late.log",
    benchmark:
        "logs/bam_correlation_pe_edu_late.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: MultiBam summary..." &&
        multiBamSummary bins \
        --bamfiles {input} \
        --minMappingQuality 20 \
        --labels late_rep1 late_uv_rep1 late_rep2 late_uv_rep2 \
        -out {output.out} --outRawCounts {output.raw_out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
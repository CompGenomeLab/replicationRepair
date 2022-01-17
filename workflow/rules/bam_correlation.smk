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
        "results/edu/R21102739-EdUrep3-UV1-5hrls3_combined/R21102739-EdUrep3-UV1-5hrls3_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21102739-EdUrep3-UV3-5hrls3_combined/R21102739-EdUrep3-UV3-5hrls3_combined_hg19_sorted_rmdup.bam",
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
        --labels early_rep1 late_rep1 early_uv_rep1 late_uv_rep1 early_rep2 late_rep2 early_uv_rep2 late_uv_rep2 early_uv_rep3 late_uv_rep3 \
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
        "results/edu/R21102739-EdUrep3-UV1-5hrls3_combined/R21102739-EdUrep3-UV1-5hrls3_combined_hg19_sorted_rmdup.bam",
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
        --labels early_rep1 early_uv_rep1 early_rep2 early_uv_rep2 early_uv_rep3 \
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
        "results/edu/R21102739-EdUrep3-UV3-5hrls3_combined/R21102739-EdUrep3-UV3-5hrls3_combined_hg19_sorted_rmdup.bam",
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
        --labels late_rep1 late_uv_rep1 late_rep2 late_uv_rep2 late_uv_rep3 \
        -out {output.out} --outRawCounts {output.raw_out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule bam_correlation_okseq:
    input:
        "results/okseq/IMR90/IMR90_se_hg19_sorted.bam",
        "results/okseq/GM06990/GM06990_se_hg19_sorted.bam",
        "results/okseq/HeLa/HeLa_se_hg19_sorted.bam",
    output:
        out="results/okseq/readCounts_okseq.npz",
        raw_out="results/okseq/readCounts_okseq.tab",
    log:
        "logs/bam_correlation_okseq.log",
    benchmark:
        "logs/bam_correlation_okseq.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: MultiBam summary..." &&
        multiBamSummary bins \
        --bamfiles {input} \
        --minMappingQuality 20 \
        --labels IMR90 GM06990 HeLa \
        -out {output.out} --outRawCounts {output.raw_out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule bam_correlation_chipseq:
    input:
        "results/chipseq/HeLa_H2AFZ_A/HeLa_H2AFZ_A_se_hg19.bam",
        "results/chipseq/HeLa_H3K27ac_A/HeLa_H3K27ac_A_se_hg19.bam",
        "results/chipseq/HeLa_H3K27me3_A/HeLa_H3K27me3_A_se_hg19.bam",
        "results/chipseq/HeLa_H3K27me3_B/HeLa_H3K27me3_B_se_hg19.bam",
        "results/chipseq/HeLa_H3K36me3_A/HeLa_H3K36me3_A_se_hg19.bam",
        "results/chipseq/HeLa_H3K36me3_B/HeLa_H3K36me3_B_se_hg19.bam",
        "results/chipseq/HeLa_H3K4me1_A/HeLa_H3K4me2_A_se_hg19.bam",
        "results/chipseq/HeLa_H3K4me3_B/HeLa_H3K4me3_B_se_hg19.bam",
        "results/chipseq/HeLa_H3K4me3_C/HeLa_H3K4me3_C_se_hg19.bam",
        "results/chipseq/HeLa_H3K79me2_A/HeLa_H3K79me2_A_se_hg19.bam",
        "results/chipseq/HeLa_H3K9ac_A/HeLa_H3K9ac_A_se_hg19.bam",
        "results/chipseq/HeLa_H3K9me3_A/HeLa_H3K9me3_A_se_hg19.bam",
        "results/chipseq/HeLa_H4K20me1_A/HeLa_H4K20me1_A_se_hg19.bam",
    output:
        out="results/chipseq/readCounts_chipseq.npz",
        raw_out="results/chipseq/readCounts_chipseq.tab",
    log:
        "logs/bam_correlation_chipseq.log",
    benchmark:
        "logs/bam_correlation_chipseq.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: MultiBam summary..." &&
        multiBamSummary bins \
        --bamfiles {input} \
        --minMappingQuality 20 \
        --labels HeLa_H2AFZ_A HeLa_H3K27ac_A HeLa_H3K27me3_A HeLa_H3K27me3_B \
            HeLa_H3K36me3_A HeLa_H3K36me3_B HeLa_H3K4me1_A HeLa_H3K4me2_A \
            HeLa_H3K4me3_B HeLa_H3K4me3_C HeLa_H3K79me2_A HeLa_H3K9ac_A \
            HeLa_H3K9me3_A HeLa_H4K20me1_A \
        -out {output.out} --outRawCounts {output.raw_out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule bam_correlation_edu_repli:
    input:
        "results/edu/R21061297-EdUrep-2hrls_combined/R21061297-EdUrep-2hrls_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21061297-EdUrep-4hrls_combined/R21061297-EdUrep-4hrls_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21061297-EdUrep-UV1-5hrls_combined/R21061297-EdUrep-UV1-5hrls_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21061297-EdUrep-UV3-5hrls_combined/R21061297-EdUrep-UV3-5hrls_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21071354-EdUrep2-2hrls2_combined/R21071354-EdUrep2-2hrls2_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21071354-EdUrep2-4hrls2_combined/R21071354-EdUrep2-4hrls2_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21071354-EdUrep2-UV1-5hrls2_combined/R21071354-EdUrep2-UV1-5hrls2_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21071354-EdUrep2-UV3-5hrls2_combined/R21071354-EdUrep2-UV3-5hrls2_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21102739-EdUrep3-UV1-5hrls3_combined/R21102739-EdUrep3-UV1-5hrls3_combined_hg19_sorted_rmdup.bam",
        "results/edu/R21102739-EdUrep3-UV3-5hrls3_combined/R21102739-EdUrep3-UV3-5hrls3_combined_hg19_sorted_rmdup.bam",
        "results/edu/repliseqG1b.bam",
        "results/edu/repliseqS1.bam",
        "results/edu/repliseqS2.bam",
        "results/edu/repliseqS3.bam",
        "results/edu/repliseqS4.bam",
        "results/edu/repliseqG2.bam",
    output:
        out="results/edu/readCounts_repli.npz",
        raw_out="results/edu/readCounts_repli.tab",
    log:
        "logs/bam_correlation_edu_repli.log",
    benchmark:
        "logs/bam_correlation_edu_repli.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: MultiBam summary..." &&
        multiBamSummary bins \
        --bamfiles {input} \
        --minMappingQuality 20 \
        --labels early_rep1 late_rep1 early_uv_rep1 late_uv_rep1 early_rep2 late_rep2 early_uv_rep2 late_uv_rep2 early_uv_rep3 late_uv_rep3 G1b S1 S2 S3 S4 G2 \
        -out {output.out} --outRawCounts {output.raw_out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
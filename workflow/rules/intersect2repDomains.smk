rule intersect2repDomains:
    input:
        rep="results/regions/repdomains_uv_mean0.5.bed",
        sns="results/regions/sns_seq_hela.bed",
        #hirfd="results/regions/hi.rfd_org.bed",
        iz="results/regions/iz_hela.bed",
        chromhmm="results/regions/wgEncodeAwgSegmentationChromhmmHelas3.bed",
        iz_comb="results/okseq/HeLa_intersect2_GM06990_IMR90.bed",
        iz_hela="results/okseq/HeLa_no_overlap.bed",
    output:
        iz="results/regions/iz_hela_repdomains_uv_mean0.5.bed",  
        sns="results/regions/sns_seq_hela_repdomains_uv_mean0.5.bed",
        #hirfd="results/regions/hi.rfd_hela_repdomains_uv_mean0.5.bed",
        chromhmm="results/regions/chromhmm_hela_repdomains_uv_mean0.5.bed",
        iz_comb_scores="results/regions/iz_hela_to_gm_imr_repdomains_uv_mean0.5_with_scores.bed",
        iz_comb="results/regions/iz_hela_to_gm_imr_repdomains_uv_mean0.5.bed",
        iz_hela="results/regions/iz_hela_no_overlap_repdomains_uv_mean0.5.bed",
    log:
        "logs/intersect2repDomains.log",
    benchmark:
        "logs/intersect2repDomains.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """  
        (echo "`date -R`: Intersecting with {input.iz}..." &&
        bedtools intersect \
        -a {input.iz} \
        -b {input.rep} -wa -wb -f 0.5 |& 
        awk '{{print $1"\\t"$2"\\t"$3"\\t""iz_"$10"\\t"$5"\\t"$6}}' \
        > {output.iz} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Intersecting with {input.sns}..." &&
        bedtools intersect \
        -a {input.sns} \
        -b {input.rep} -wa -wb -f 0.5 |& 
        awk '{{print $1"\\t"$2"\\t"$3"\\t""sns_"$10"\\t"$5"\\t"$6}}' \
        > {output.sns} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Intersecting with {input.chromhmm}..." &&
        bedtools intersect \
        -a {input.chromhmm} \
        -b {input.rep} -wa -wb -f 0.5 |& 
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"".""\\t"$13}}' \
        > {output.chromhmm} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Intersecting with {input.iz_comb}..." &&
        bedtools intersect \
        -a {input.iz_comb} \
        -b {input.rep} -wa -wb -f 0.5 |& 
        awk '{{print $1"\\t"$2"\\t"$3"\\t""HeLa_"$9"\\t"".""\\t""."}}' \
        > {output.iz_comb} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Intersecting with {input.iz_comb} with scores..." &&
        bedtools intersect \
        -a {input.iz_comb} \
        -b {input.rep} -wa -wb -f 0.5 |& 
        awk '{{print $1"\\t"$2"\\t"$3"\\t""HeLa_"$9"_"$5"\\t"".""\\t""."}}' \
        > {output.iz_comb_scores} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Intersecting with {input.iz_hela}..." &&
        bedtools intersect \
        -a {input.iz_hela} \
        -b {input.rep} -wa -wb -f 0.5 |& 
        awk '{{print $1"\\t"$2"\\t"$3"\\t""HeLa_"$9"\\t"".""\\t""."}}' \
        > {output.iz_hela} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule intersect2repDomains_methyl:
    input:
        rep="results/regions/repdomains_uv_mean0.5.bed",
        methyl="results/regions/methylation_shuf_1m_windows_1_3.bed",
    output:
        "results/regions/methylation_shuf_1m_repdomains_windows_1_3.bed",
    log:
        "logs/intersect2repDomains_methyl.log",
    benchmark:
        "logs/intersect2repDomains_methyl.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """  
        (echo "`date -R`: Intersecting with {input.methyl}..." &&
        bedtools intersect \
        -a {input.methyl} \
        -b {input.rep} -wa -wb -f 0.5 |& 
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"_"$10"\\t"$5"\\t"$6}}' \
        > {output} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

rule intersect2repDomains_hirfd:
    input:
        rep="results/regions/repdomains_uv_mean0.5.bed",
        hirfd="results/regions/hi.rfd_org.bed",
    output:
        hirfd="results/regions/hi.rfd_hela_repdomains_uv_mean0.5.bed",
    log:
        "logs/intersect2repDomains_hirfd.log",
    benchmark:
        "logs/intersect2repDomains_hirfd.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """  
        (echo "`date -R`: Intersecting with {input.hirfd}..." &&
        bedtools intersect \
        -a {input.hirfd} \
        -b {input.rep} -wa -wb -f 0.5 |& 
        awk '{{print $1"\\t"$2"\\t"$3"\\t""HighRFD_"$10"\\t"$5"\\t"$6}}' \
        > {output.hirfd} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """
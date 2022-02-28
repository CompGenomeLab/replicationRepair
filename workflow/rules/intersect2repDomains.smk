rule intersect2repDomains:
    input:
        rep="results/regions/repdomains_uv_mean0.5.bed",
        iz="resources/samples/iz_hela.bed",
        chromhmm="resources/samples/chromHMM/wgEncodeAwgSegmentationChromhmmHelas3.bed",
        iz_comb="results/okseq/HeLa_intersect2_GM06990_IMR90.bed",
        iz_hela="results/okseq/HeLa_no_overlap.bed",
    output:
        iz="results/regions/iz_hela_repdomains_uv_mean0.5.bed",  
        chromhmm="results/regions/chromhmm_hela_repdomains_uv_mean0.5.bed",
        iz_comb_scores="results/regions/iz_hela_to_gm_imr_repdomains_uv_mean0.5_with_scores.bed",
        iz_comb="results/regions/iz_hela_to_gm_imr_repdomains_uv_mean0.5.bed",
        iz_hela="results/regions/iz_hela_no_overlap_repdomains_uv_mean0.5.bed",
    log:
        "logs/rule/analysis/intersect2repDomains.log",
    benchmark:
        "logs/rule/analysis/intersect2repDomains.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
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
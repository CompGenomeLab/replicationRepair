rule intersect2repDomains:
    input:
        rep="results/regions/repdomains_uv_mean0.5.bed",
        sns="results/regions/sns_seq_hela.bed",
        iz="results/regions/iz_hela.bed",
        chromhmm="results/regions/wgEncodeAwgSegmentationChromhmmHelas3.bed"
    output:
        iz="results/regions/iz_hela_repdomains_uv_mean0.5.bed",  
        sns="results/regions/sns_seq_hela_repdomains_uv_mean0.5.bed",
        chromhmm="results/regions/chromhmm_hela_repdomains_uv_mean0.5.bed"
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
        """
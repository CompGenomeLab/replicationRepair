rule okseqIntersect:
    input:
        hela="results/okseq/HeLa/HeLa_hg19_HMMsegments_IZ.bed",
        gm06990="results/okseq/GM06990/GM06990_hg19_HMMsegments_IZ.bed",
        imr90="results/okseq/IMR90/IMR90_hg19_HMMsegments_IZ.bed",
    output:
        hela="results/okseq/HeLa_intersect2_GM06990_IMR90.txt",
        gm06990="results/okseq/GM06990_intersect2_HeLa_IMR90.txt",
        imr90="results/okseq/IMR90_intersect2_GM06990_HeLa.txt",
    log:
        "logs/okseqIntersect.log",
    benchmark:
        "logs/okseqIntersect.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """     
        (echo "`date -R`: Intersecting to HeLa..." &&
        bedtools intersect -wao \
        -a {input.hela} \
        -b {input.gm06990} {input.imr90} \
        -names GM06990 IMR90 \
        > {output.hela} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Intersecting to GM06990..." &&
        bedtools intersect -wao \
        -a {input.gm06990} \
        -b {input.hela} {input.imr90} \
        -names HeLa IMR90 \
        > {output.gm06990} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Intersecting to IMR90..." &&
        bedtools intersect -wao \
        -a {input.imr90} \
        -b {input.gm06990} {input.hela} \
        -names GM06990 HeLa \
        > {output.imr90} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1        
        """
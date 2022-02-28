rule ts_nts_xr:
    input:
        plus="resources/samples/XR/{samples}_{build}_sorted_plus.bed",
        minus="resources/samples/XR/{samples}_{build}_sorted_minus.bed",
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes_knownCanonical_stranded.bed",
    output:
        comb=temp("results/XR/{samples}/{samples}_{build}_sorted_chr.bed"),
        TS=temp("results/XR/{samples}/{samples}_{build}_sorted_TS.bed"),
        NTS=temp("results/XR/{samples}/{samples}_{build}_sorted_NTS.bed"),
        TSNTS="results/XR/{samples}/{samples}_{build}_sorted_TSNTS.bed",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_ts_nts_xr.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_ts_nts_xr.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:
        """
        (echo "`date -R`: Combine strands..." &&
        cat {input.plus} {input.minus} | \
        sort -k1,1 -k2,2n -k3,3n > {output.comb} &&
        echo "`date -R`: Success! Strands combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Counting reads on TS..." &&
        bedtools intersect \
        -sorted -a {input.genes} \
        -b {output.comb} \
        -wa -c -S -F 0.5 > {output.TS} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Counting reads on NTS..." &&
        bedtools intersect \
        -sorted -a {input.genes} \
        -b {output.comb} \
        -wa -c -s -F 0.5 > {output.NTS} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Combine TS & NTS counts..." &&
        paste {output.TS} {output.NTS} | \
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$6"\\t"$7"\\t"$14}}' > \
        {output.TSNTS} &&
        echo "`date -R`: Success! TS & NTS counts are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule ts_nts_ds:
    input:
        plus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes_knownCanonical_stranded.bed",
    output:
        comb=temp("results/DS/{samples}/{samples}_{build}_sorted_dipyrimidines_chr.bed"),
        TS=temp("results/DS/{samples}/{samples}_{build}_sorted_dipyrimidines_TS.bed"),
        NTS=temp("results/DS/{samples}/{samples}_{build}_sorted_dipyrimidines_NTS.bed"),
        TSNTS="results/DS/{samples}/{samples}_{build}_sorted_dipyrimidines_TSNTS.bed",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_ts_nts_ds.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_ts_nts_ds.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:
        """
        (echo "`date -R`: Combine strands..." &&
        cat {input.plus} {input.minus} | \
        sort -k1,1 -k2,2n -k3,3n > {output.comb} &&
        echo "`date -R`: Success! Strands combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Counting reads on TS..." &&
        bedtools intersect \
        -sorted -a {input.genes} \
        -b {output.comb} \
        -wa -c -S -F 0.5 > {output.TS} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Counting reads on NTS..." &&
        bedtools intersect \
        -sorted -a {input.genes} \
        -b {output.comb} \
        -wa -c -s -F 0.5 > {output.NTS} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Combine TS & NTS counts..." &&
        paste {output.TS} {output.NTS} | \
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$6"\\t"$7"\\t"$14}}' > \
        {output.TSNTS} &&
        echo "`date -R`: Success! TS & NTS counts are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
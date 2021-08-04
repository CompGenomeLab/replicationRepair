rule countMotifs:
    input:
        bed="results/regions/{regions}.bed",
        genome="resources/ref_genomes/hg19/genome_hg19.fa",
    output:
        fa=temp("results/regions/{regions}.fa"),
        counts="results/mutation/{regions}_counts.txt",
    params:
        "tc cc ga gg",
    log:
        "logs/mutation/countMotifs_{regions}.log",
    benchmark:
        "logs/mutation/countMotifs_{regions}.benchmark.txt",
    conda:
        "../envs/countMotifs.yaml"
    shell:
        """
        (echo "`date -R`: Converting bed region file to fasta..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output.fa} \
        -name &&
        echo "`date -R`: Success!" ||
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Counting the given motif..." &&
        workflow/scripts/countMotif.py \
        -i {output.fa} \
        -o {output.counts} \
        -m {params} \
        --agg &&
        echo "`date -R`: Success! Counting is done." ||
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

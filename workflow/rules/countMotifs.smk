rule countMotifs:
    input:
        bed="results/regions/{regions}.bed",
        genome="resources/ref_genomes/hg19/genome_hg19.fa",
    output:
        fa=temp("results/regions/{regions}.fa"),
        counts="results/regions/{regions}_counts.txt",
    params:
        "tc cc ga gg",
    log:
        "logs/regions/countMotifs_{regions}.log",
    benchmark:
        "logs/regions/countMotifs_{regions}.benchmark.txt",
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
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Counting the given motif..." &&
        python3 workflow/scripts/countMotif.py \
        -i {output.fa} \
        -o {output.counts} \
        -m {params} \
        --agg &&
        echo "`date -R`: Success! Counting is done." ||
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule countMotifs_kmer:
    input:
        bed="results/regions/{regions}.bed",
        genome="resources/ref_genomes/hg19/genome_hg19.fa",
    output:
        plus=temp("results/regions/{regions}_plus.bed"),
        minus=temp("results/regions/{regions}_minus.bed"),
        stranded="results/regions/{regions}_stranded.bed",
        fa=temp("results/regions/{regions}_stranded.fa"),
        counts="results/regions/{regions}_kmer_counts.txt",
    params:
        kmer=5,
    log:
        "logs/regions/countMotifs_kmer_{regions}.log",
    benchmark:
        "logs/regions/countMotifs_kmer_{regions}.benchmark.txt",
    conda:
        "../envs/countMotifs.yaml"
    shell:
        """
        (echo "`date -R`: Make bed file stranded..." &&
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"_+""\\t"$5"\\t""+"}}' {input.bed} \
        > {output.plus} &&
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"_-""\\t"$5"\\t""-"}}' {input.bed} \
        > {output.minus} &&
        cat {output.plus} {output.minus} > {output.stranded} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Converting bed region file to fasta..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {output.stranded} \
        -fo {output.fa} \
        -name \
        -s &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Counting until the {params.kmer}-mer..." &&
        python3 workflow/scripts/countByKmer.py \
        -i {output.fa} \
        -o {output.counts} \
        --kmer "{params.kmer}" \
        --agg &&
        echo "`date -R`: Success! Counting is done." ||
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule countMotifs_genome:
    input:
        genome="resources/ref_genomes/hg19/genome_hg19.fa",
    output:
        genome="results/regions/genome_T_counts.txt",
    log:
        "logs/regions/countMotifs_genome.log",
    benchmark:
        "logs/regions/countMotifs_genome.benchmark.txt",
    conda:
        "../envs/countMotifs.yaml"
    shell:
        """
        (echo "`date -R`: Counting the given motif..." &&
        python3 workflow/scripts/countMotif.py \
        -i {input.genome} \
        -o {output.genome} \
        -m "T" \
        --agg &&
        echo "`date -R`: Success! Counting is done." ||
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """
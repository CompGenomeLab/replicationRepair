rule genomecov_edu:
    input:
        plus="results/edu/{samples}/{samples}_{build}_sorted_plus.bed",
        minus="results/edu/{samples}/{samples}_{build}_sorted_minus.bed",
        ref_genome="resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        plus=temp("results/edu/{samples}/{samples}_{build}_sorted_plus.bdg"),
        minus=temp("results/edu/{samples}/{samples}_{build}_sorted_minus.bdg"),
    params:
        read=lambda w, input: mappedReads(input[0], input[1]),
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_genomecov_edu.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_genomecov_edu.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Calculating genome coverage of {input.plus}..." &&
        bedtools genomecov \
        -i {input.plus} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.plus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1


        (echo "`date -R`: Calculating genome coverage of {input.minus}..." &&
        bedtools genomecov \
        -i {input.minus} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print -1000000/$1}}') \
        > {output.minus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule genomecov_edu_v2:
    input:
        bed="results/edu/{samples}/{samples}_{build}_pe.bed",
        genome="resources/ref_genomes/hg19/genome_hg19_50kb.bed"
    output:
        "results/edu/{samples}/{samples}_{build}.bdg",
    params:
        read=lambda w, input: mappedReads(input[0]),
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_genomecov_edu_v2.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_genomecov_edu_v2.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Calculating genome coverage..." &&
        bedtools intersect -c \
        -b {input.bed} \
        -a {input.genome} |&
        awk '{{print $1,$2,$3,$4*1e+06/{params.read} }}' OFS='\\t' | grep chr > {output} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """
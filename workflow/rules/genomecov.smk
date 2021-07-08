
rule genomecov_input:
    input:
        plus="results/input/{samples}/{samples}_{build}_sorted_plus.bed",
        minus="results/input/{samples}/{samples}_{build}_sorted_minus.bed",
        ref_genome="resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        plus=temp("results/input/{samples}/{samples}_{build}_sorted_plus.bdg"),
        minus=temp("results/input/{samples}/{samples}_{build}_sorted_minus.bdg"),
    params:
        read=lambda w, input: mappedReads(input),
    log:
        "logs/{samples}/{samples}_{build}_genomecov_input.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_genomecov_input.benchmark.txt",
    conda:
        "../envs/genomecov.yaml"
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
        echo "`date -R`: Process failed...") > {log} 2>&1


        (echo "`date -R`: Calculating genome coverage of {input.minus}..." &&
        bedtools genomecov \
        -i {input.minus} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print -1000000/$1}}') \
        > {output.minus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule genomecov_edu:
    input:
        plus="results/edu/{samples}/{samples}_{build}_sorted_plus.bed",
        minus="results/edu/{samples}/{samples}_{build}_sorted_minus.bed",
        ref_genome="resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        plus=temp("results/edu/{samples}/{samples}_{build}_sorted_plus.bdg"),
        minus=temp("results/edu/{samples}/{samples}_{build}_sorted_minus.bdg"),
    params:
        read=lambda w, input: mappedReads(input),
    log:
        "logs/{samples}/{samples}_{build}_genomecov_edu.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_genomecov_edu.benchmark.txt",
    conda:
        "../envs/genomecov.yaml"
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
        echo "`date -R`: Process failed...") > {log} 2>&1


        (echo "`date -R`: Calculating genome coverage of {input.minus}..." &&
        bedtools genomecov \
        -i {input.minus} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print -1000000/$1}}') \
        > {output.minus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """
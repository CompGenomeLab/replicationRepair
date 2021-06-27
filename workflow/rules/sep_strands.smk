
rule sep_strands_input:
    input:
        "results/input/{samples}/{samples}_{build}_sorted_chr.bed",
    output:
        plus="results/input/{samples}/{samples}_{build}_sorted_plus.bed",
        minus="results/input/{samples}/{samples}_{build}_sorted_minus.bed",
    log:
        "logs/{samples}/{samples}_{build}_sep_strands.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_sep_strands.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Separating plus stranded reads..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus} &&
        echo "`date -R`: Success! Reads are separated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Separating minus stranded reads..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Reads are separated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """


rule sep_strands_mutation:
    input:
        "results/mutation/{samples}/{samples}_target_mut.tsv",
    output:
        plus="results/mutation/{samples}/{samples}_target_mut_plus.tsv",
        minus="results/mutation/{samples}/{samples}_target_mut_minus.tsv",
    log:
        "logs/{samples}/{samples}_sep_strands.log",
    benchmark:
        "logs/{samples}/{samples}_sep_strands.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Separating plus stranded mutations..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus} &&
        echo "`date -R`: Success! Mutations are separated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Separating minus stranded mutations..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Mutations are separated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """


rule sep_strands_sim:
    input:
        "resources/samples/sim/{samples}_{build}_{method}_sim.bed",
    output:
        plus=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_plus.bed"),
        minus=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_minus.bed"),
    log:
        "logs/{samples}/{samples}_{build}_{method}_sim_sep_strands.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_{method}_sim_sep_strands.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Separating plus stranded mutations..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus} &&
        echo "`date -R`: Success! Mutations are separated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Separating minus stranded mutations..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Mutations are separated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """
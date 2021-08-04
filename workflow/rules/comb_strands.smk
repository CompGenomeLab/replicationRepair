rule comb_strands_mutation:
    input:
        plus="results/mutation/{samples}/{samples}_target_mut_plus_{regions}.txt",
        minus="results/mutation/{samples}/{samples}_target_mut_plus_{regions}.txt",
    output:
        comb=temp("results/mutation/{samples}/{samples}_target_mut_comb_{regions}.txt"),
        org="results/mutation/{samples}/{samples}_target_mut_comb_{regions}_org.txt",
    log:
        "logs/{samples}/{samples}_{regions}_comb_strands_mutation.log",
    benchmark:
        "logs/{samples}/{samples}_{regions}_comb_strands_mutation.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Combining plus stranded mutations..." &&
        awk '{{print $0"\\t""+"}}' {input.plus} > {output.comb} &&
        echo "`date -R`: Success! Mutations are combined." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Combining minus stranded mutations..." &&
        awk '{{print $0"\\t""-"}}' {input.minus} >> {output.comb} &&
        echo "`date -R`: Success! Mutations are combined." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Organizing..." &&
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$6"\\t"$8"\\t"$7}}' \
        {output.comb} > {output.org} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule comb_strands_mutation_intergenic:
    input:
        plus="results/mutation/{samples}/{samples}_target_mut_plus_{regions}_intergenic.txt",
        minus="results/mutation/{samples}/{samples}_target_mut_minus_{regions}_intergenic.txt",
    output:
        comb=temp("results/mutation/{samples}/{samples}_target_mut_comb_{regions}_intergenic.txt"),
        org="results/mutation/{samples}/{samples}_target_mut_comb_{regions}_intergenic_org.txt",
    log:
        "logs/{samples}/{samples}_{regions}_comb_strands_mutation_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{regions}_comb_strands_mutation_intergenic.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Combining plus stranded mutations..." &&
        awk '{{print $0"\\t""+"}}' {input.plus} > {output.comb} &&
        echo "`date -R`: Success! Mutations are combined." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Combining minus stranded mutations..." &&
        awk '{{print $0"\\t""-"}}' {input.minus} >> {output.comb} &&
        echo "`date -R`: Success! Mutations are combined." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Organizing..." &&
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$6"\\t"$8"\\t"$7}}' \
        {output.comb} > {output.org} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """
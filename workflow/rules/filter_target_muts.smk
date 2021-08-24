
rule filter_target_muts:
    input:
        "results/mutation/{samples}/{samples}.bed",
    output:
        "results/mutation/{samples}/{samples}_target_mut.tsv",   
    log:
        "logs/{samples}/{samples}_filter_target_muts.log",
    benchmark:
        "logs/{samples}/{samples}_filter_target_muts.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Filtering target mutations..." &&
        grep -e 'C_T_[TC]' -e 'G_A_..[GA]' {input} | \
        awk -v t="\\t" '{{\
        if ($4 ~ "C_T"){{print $0t"0"t"+"}}; \
        if ($4 ~ "G_A"){{print $0t"0"t"-"}}}}' > \
        {output} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

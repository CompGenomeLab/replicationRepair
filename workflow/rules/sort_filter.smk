
rule sort_filter_input:
    input:
        lambda w: input4filter(w, config["sample_input"], "input"),
    output:
        "results/input/{samples}/{samples}_{build}_sorted_chr.bed",
    params:
        filt="'^chr([1-9]|1[0-9]|2[0-2]|X)'",
    log:
        "logs/{samples}/{samples}_{build}_sort_filter.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_sort_filter.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Sorting and filtering bed file by chromosomes..." &&
        sort -u -k1,1 -k2,2n -k3,3n {input} |&
        egrep {params.filt} > {output} &&
        echo "`date -R`: Success! Bed file is filtered." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """
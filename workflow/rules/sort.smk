rule sort_xr:
    input:
        plus="resources/samples/XR/{samples}_{build}_sorted_plus.bed",
        minus="resources/samples/XR/{samples}_{build}_sorted_minus.bed",
    output:
        plus=temp("results/XR/{samples}/{samples}_{build}_xr_plus_sorted.txt"),
        minus=temp("results/XR/{samples}/{samples}_{build}_xr_minus_sorted.txt"),
    params:
        filt="'^chr([1-9]|1[0-9]|2[0-2]|X)'",
    log:
        "logs/{samples}/{samples}_{build}_sort_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_sort_xr.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Sorting and filtering bed file by chromosomes..." &&
        sort -u -k1,1 -k2,2n -k3,3n {input.plus} |&
        egrep {params.filt} > {output.plus} &&
        echo "`date -R`: Success! Bed file is filtered." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Sorting and filtering bed file by chromosomes..." &&
        sort -u -k1,1 -k2,2n -k3,3n {input.minus} |&
        egrep {params.filt} > {output.minus} &&
        echo "`date -R`: Success! Bed file is filtered." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule sort_ds:
    input:
        plus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
    output:
        plus=temp("results/DS/{samples}/{samples}_{build}_ds_dipyrimidines_plus_sorted.txt"),
        minus=temp("results/DS/{samples}/{samples}_{build}_ds_dipyrimidines_minus_sorted.txt"),   
    params:
        filt="'^chr([1-9]|1[0-9]|2[0-2]|X)'",
    log:
        "logs/{samples}/{samples}_{build}_sort_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_sort_ds.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Sorting and filtering bed file by chromosomes..." &&
        sort -u -k1,1 -k2,2n -k3,3n {input.plus} |&
        egrep {params.filt} > {output.plus} &&
        echo "`date -R`: Success! Bed file is filtered." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Sorting and filtering bed file by chromosomes..." &&
        sort -u -k1,1 -k2,2n -k3,3n {input.minus} |&
        egrep {params.filt} > {output.minus} &&
        echo "`date -R`: Success! Bed file is filtered." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule sort_sim:
    input:
        plus_sim="results/sim/{samples}/{samples}_{build}_{method}_sim_plus.bed",
        minus_sim="results/sim/{samples}/{samples}_{build}_{method}_sim_minus.bed",
    output:
        plus=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_plus_sorted.txt"),
        minus=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_minus_sorted.txt"),
    params:
        filt="'^chr([1-9]|1[0-9]|2[0-2]|X)'",
    log:
        "logs/{samples}/{samples}_{build}_sort_sim_{method}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_sort_sim_{method}.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Sorting and filtering bed file by chromosomes..." &&
        sort -u -k1,1 -k2,2n -k3,3n {input.plus} |&
        egrep {params.filt} > {output.plus} &&
        echo "`date -R`: Success! Bed file is filtered." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Sorting and filtering bed file by chromosomes..." &&
        sort -u -k1,1 -k2,2n -k3,3n {input.minus} |&
        egrep {params.filt} > {output.minus} &&
        echo "`date -R`: Success! Bed file is filtered." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """
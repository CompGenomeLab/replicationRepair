rule more_info_xr:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}_combined.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}_combined.txt",
        reads_p="resources/samples/XR/{samples}_{build}_sorted_plus.bed",
        reads_m="resources/samples/XR/{samples}_{build}_sorted_minus.bed",
    output:
        plus=temp("results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}_combined_info.txt"),
        minus=temp("results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}_combined_info.txt"),
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[2], input[3]),
    log:
        "logs/{samples}/{samples}_{build}_more_info_{regions}_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_more_info_{regions}_xr.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Add info (plus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.plus} \
        -o {output.plus} \
        -c {params.info} "." {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} "." {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule more_info_ds:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}_combined.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}_combined.txt",
        reads_p="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        reads_m="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
    output:
        plus=temp("results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}_combined_info.txt"),
        minus=temp("results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}_combined_info.txt"),
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[2], input[3]),
    log:
        "logs/{samples}/{samples}_{build}_more_info_{regions}_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_more_info_{regions}_ds.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Add info (plus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.plus} \
        -o {output.plus} \
        -c {params.info} "." {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} "." {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """
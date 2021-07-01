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
        -c {params.info} + {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} - {params.mappedReads} &&
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
        -c {params.info} + {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} - {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule more_info_sim:
    input:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{regions}_combined.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{regions}_combined.txt",
        reads_p="results/sim/{samples}/{samples}_{build}_{method}_sim_plus.bed",
        reads_m="results/sim/{samples}/{samples}_{build}_{method}_sim_minus.bed",
    output:
        plus=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{regions}_combined_info.txt"),
        minus=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{regions}_combined_info.txt"),
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[2], input[3]),
    log:
        "logs/{samples}/{samples}_{build}_more_info_{method}_sim_{regions}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_more_info_{method}_sim_{regions}.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Add info (plus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.plus} \
        -o {output.plus} \
        -c {params.info} + {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} - {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule more_info_xr_intergenic:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}_intergenic_combined.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}_intergenic_combined.txt",
        reads_p="results/XR/{samples}/{samples}_{build}_intergenic_sorted_plus.bed",
        reads_m="results/XR/{samples}/{samples}_{build}_intergenic_sorted_minus.bed",
    output:
        plus=temp("results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}_intergenic_combined_info.txt"),
        minus=temp("results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}_intergenic_combined_info.txt"),
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[2], input[3]),
    log:
        "logs/{samples}/{samples}_{build}_more_info_{regions}_xr_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_more_info_{regions}_xr_intergenic.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Add info (plus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.plus} \
        -o {output.plus} \
        -c {params.info} + {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} - {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule more_info_ds_intergenic:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}_intergenic_combined.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}_intergenic_combined.txt",
        reads_p="results/DS/{samples}/{samples}_{build}_intergenic_sorted_ds_dipyrimidines_plus.bed",
        reads_m="results/DS/{samples}/{samples}_{build}_intergenic_sorted_ds_dipyrimidines_minus.bed",
    output:
        plus=temp("results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}_intergenic_combined_info.txt"),
        minus=temp("results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}_intergenic_combined_info.txt"),
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[2], input[3]),
    log:
        "logs/{samples}/{samples}_{build}_more_info_{regions}_ds_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_more_info_{regions}_ds_intergenic.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Add info (plus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.plus} \
        -o {output.plus} \
        -c {params.info} + {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} - {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule more_info_sim_intergenic:
    input:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{regions}_intergenic_combined.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{regions}_intergenic_combined.txt",
        reads_p="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_intergenic.bed",
        reads_m="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_intergenic.bed",
    output:
        plus=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{regions}_intergenic_combined_info.txt"),
        minus=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{regions}_intergenic_combined_info.txt"),
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[2], input[3]),
    log:
        "logs/{samples}/{samples}_{build}_more_info_{method}_sim_{regions}_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_more_info_{method}_sim_{regions}_intergenic.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Add info (plus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.plus} \
        -o {output.plus} \
        -c {params.info} + {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} - {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule more_info_noWindows_xr:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{noWindows}.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{noWindows}.txt",
        reads_p="resources/samples/XR/{samples}_{build}_sorted_plus.bed",
        reads_m="resources/samples/XR/{samples}_{build}_sorted_minus.bed",
    output:
        plus=temp("results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{noWindows}_info.txt"),
        minus=temp("results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{noWindows}_info.txt"),
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[2], input[3]),
    log:
        "logs/{samples}/{samples}_{build}_more_info_noWindows_{noWindows}_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_more_info_noWindows_{noWindows}_xr.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Add info (plus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.plus} \
        -o {output.plus} \
        -c {params.info} + {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} - {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule more_info_noWindows_ds:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{noWindows}.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{noWindows}.txt",
        reads_p="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        reads_m="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
    output:
        plus=temp("results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{noWindows}_info.txt"),
        minus=temp("results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{noWindows}_info.txt"),
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[2], input[3]),
    log:
        "logs/{samples}/{samples}_{build}_more_info_noWindows_{noWindows}_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_more_info_noWindows_{noWindows}_ds.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Add info (plus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.plus} \
        -o {output.plus} \
        -c {params.info} + {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} - {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule more_info_noWindows_sim:
    input:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{noWindows}.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{noWindows}.txt",
        reads_p="results/sim/{samples}/{samples}_{build}_{method}_sim_plus.bed",
        reads_m="results/sim/{samples}/{samples}_{build}_{method}_sim_minus.bed",
    output:
        plus=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{noWindows}_info.txt"),
        minus=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{noWindows}_info.txt"),
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[2], input[3]),
    log:
        "logs/{samples}/{samples}_{build}_more_info_noWindows_{method}_sim_{noWindows}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_more_info_noWindows_{method}_sim_{noWindows}.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Add info (plus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.plus} \
        -o {output.plus} \
        -c {params.info} + {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} - {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule more_info_noWindows_xr_intergenic:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{noWindows}_intergenic.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{noWindows}_intergenic.txt",
        reads_p="results/XR/{samples}/{samples}_{build}_intergenic_sorted_plus.bed",
        reads_m="results/XR/{samples}/{samples}_{build}_intergenic_sorted_minus.bed",
    output:
        plus=temp("results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{noWindows}_intergenic_info.txt"),
        minus=temp("results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{noWindows}_intergenic_info.txt"),
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[2], input[3]),
    log:
        "logs/{samples}/{samples}_{build}_more_info_noWindows_{noWindows}_xr_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_more_info_noWindows_{noWindows}_xr_intergenic.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Add info (plus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.plus} \
        -o {output.plus} \
        -c {params.info} + {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} - {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule more_info_noWindows_ds_intergenic:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{noWindows}_intergenic.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{noWindows}_intergenic.txt",
        reads_p="results/DS/{samples}/{samples}_{build}_intergenic_sorted_ds_dipyrimidines_plus.bed",
        reads_m="results/DS/{samples}/{samples}_{build}_intergenic_sorted_ds_dipyrimidines_minus.bed",
    output:
        plus=temp("results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{noWindows}_intergenic_info.txt"),
        minus=temp("results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{noWindows}_intergenic_info.txt"),
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[2], input[3]),
    log:
        "logs/{samples}/{samples}_{build}_more_info_noWindows_{noWindows}_ds_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_more_info_noWindows_{noWindows}_ds_intergenic.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Add info (plus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.plus} \
        -o {output.plus} \
        -c {params.info} + {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} - {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule more_info_noWindows_sim_intergenic:
    input:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{noWindows}_intergenic.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{noWindows}_intergenic.txt",
        reads_p="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_intergenic.bed",
        reads_m="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_intergenic.bed",
    output:
        plus=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{noWindows}_intergenic_info.txt"),
        minus=temp("results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{noWindows}_intergenic_info.txt"),
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[2], input[3]),
    log:
        "logs/{samples}/{samples}_{build}_more_info_noWindows_{method}_sim_{noWindows}_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_more_info_noWindows_{method}_sim_{noWindows}_intergenic.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Add info (plus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.plus} \
        -o {output.plus} \
        -c {params.info} + {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Add info (minus strand)..." &&
        python3 workflow/scripts/addColumns.py \
        -i {input.minus} \
        -o {output.minus} \
        -c {params.info} - {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

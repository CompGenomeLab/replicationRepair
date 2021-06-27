rule combine_windows_xr:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}.txt",
    output:
        plus=temp("results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}_combined.txt"),
        minus=temp("results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}_combined.txt"),
    params:
        lambda w: getCombine(w.regions, config["region_comb_opt"], config["regions"]), 
    log:
        "logs/{samples}/{samples}_{build}_combine_windows_{regions}_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_combine_windows_{regions}_xr.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combine windows (plus strand)..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {input.plus} \
        {params} \
        -o {output.plus} &&
        echo "`date -R`: Success! Windows are combined." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Combine windows (minus strand)..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {input.minus} \
        {params} \
        -o {output.minus} &&
        echo "`date -R`: Success! Windows are combined." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule combine_windows_ds:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}.txt",
    output:
        plus=temp("results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}_combined.txt"),
        minus=temp("results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}_combined.txt"),
    params:
        lambda w: getCombine(w.regions, config["region_comb_opt"], config["regions"]), 
    log:
        "logs/{samples}/{samples}_{build}_combine_windows_{regions}_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_combine_windows_{regions}_ds.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combine windows (plus strand)..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {input.plus} \
        {params} \
        -o {output.plus} &&
        echo "`date -R`: Success! Windows are combined." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Combine windows (minus strand)..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {input.minus} \
        {params} \
        -o {output.minus} &&
        echo "`date -R`: Success! Windows are combined." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule combine_windows_xr_sim:
    input:
        plus="results/XR/{samples}/{samples}_{build}_xr_sim_plus_{regions}.txt",
        minus="results/XR/{samples}/{samples}_{build}_xr_sim_minus_{regions}.txt",
    output:
        plus=temp("results/XR/{samples}/{samples}_{build}_xr_sim_plus_{regions}_combined.txt"),
        minus=temp("results/XR/{samples}/{samples}_{build}_xr_sim_minus_{regions}_combined.txt"),
    params:
        lambda w: getCombine(w.regions, config["region_comb_opt"], config["regions"]), 
    log:
        "logs/{samples}/{samples}_{build}_combine_windows_xr_sim_{regions}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_combine_windows_xr_sim_{regions}.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combine windows (plus strand)..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {input.plus} \
        {params} \
        -o {output.plus} &&
        echo "`date -R`: Success! Windows are combined." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Combine windows (minus strand)..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {input.minus} \
        {params} \
        -o {output.minus} &&
        echo "`date -R`: Success! Windows are combined." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule combine_windows_ds_sim:
    input:
        plus="results/DS/{samples}/{samples}_{build}_ds_sim_plus_{regions}.txt",
        minus="results/DS/{samples}/{samples}_{build}_ds_sim_minus_{regions}.txt",
    output:
        plus=temp("results/DS/{samples}/{samples}_{build}_ds_sim_plus_{regions}_combined.txt"),
        minus=temp("results/DS/{samples}/{samples}_{build}_ds_sim_minus_{regions}_combined.txt"),
    params:
        lambda w: getCombine(w.regions, config["region_comb_opt"], config["regions"]), 
    log:
        "logs/{samples}/{samples}_{build}_combine_windows_{regions}_ds_sim.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_combine_windows_{regions}_ds_sim.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combine windows (plus strand)..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {input.plus} \
        {params} \
        -o {output.plus} &&
        echo "`date -R`: Success! Windows are combined." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Combine windows (minus strand)..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {input.minus} \
        {params} \
        -o {output.minus} &&
        echo "`date -R`: Success! Windows are combined." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """
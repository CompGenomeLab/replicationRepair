rule rpkm_xr:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}_combined_info.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}_combined_info.txt",
    output:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}_combined_rpkm.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}_combined_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_{regions}_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_{regions}_xr.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Calculating RPKM values (plus strand)..." &&
        python3 workflow/scripts/RPKM.py \
        -i {input.plus} \
        -chse 2 3 \
        -c 7 \
        -mr 0 \
        -o {output.plus} &&
        echo "`date -R`: Success! RPKMs are calculated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Calculating RPKM values (minus strand)..." &&
        python3 workflow/scripts/RPKM.py \
        -i {input.minus} \
        -chse 2 3 \
        -c 7 \
        -mr 0 \
        -o {output.minus} &&
        echo "`date -R`: Success! RPKMs are calculated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule rpkm_ds:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}_combined_info.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}_combined_info.txt",
    output:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}_combined_rpkm.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}_combined_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_{regions}_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_{regions}_ds.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Calculating RPKM values (plus strand)..." &&
        python3 workflow/scripts/RPKM.py \
        -i {input.plus} \
        -chse 2 3 \
        -c 7 \
        -mr 0 \
        -o {output.plus} &&
        echo "`date -R`: Success! RPKMs are calculated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Calculating RPKM values (minus strand)..." &&
        python3 workflow/scripts/RPKM.py \
        -i {input.minus} \
        -chse 2 3 \
        -c 7 \
        -mr 0 \
        -o {output.minus} &&
        echo "`date -R`: Success! RPKMs are calculated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule rpkm_xr_sim:
    input:
        plus="results/XR/{samples}/{samples}_{build}_xr_sim_plus_{regions}_combined_info.txt",
        minus="results/XR/{samples}/{samples}_{build}_xr_sim_minus_{regions}_combined_info.txt",
    output:
        plus="results/XR/{samples}/{samples}_{build}_xr_sim_plus_{regions}_combined_rpkm.txt",
        minus="results/XR/{samples}/{samples}_{build}_xr_sim_minus_{regions}_combined_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_xr_sim_{regions}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_xr_sim_{regions}.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Calculating RPKM values (plus strand)..." &&
        python3 workflow/scripts/RPKM.py \
        -i {input.plus} \
        -chse 2 3 \
        -c 7 \
        -mr 0 \
        -o {output.plus} &&
        echo "`date -R`: Success! RPKMs are calculated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Calculating RPKM values (minus strand)..." &&
        python3 workflow/scripts/RPKM.py \
        -i {input.minus} \
        -chse 2 3 \
        -c 7 \
        -mr 0 \
        -o {output.minus} &&
        echo "`date -R`: Success! RPKMs are calculated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule rpkm_ds_sim:
    input:
        plus="results/DS/{samples}/{samples}_{build}_ds_sim_plus_{regions}_combined_info.txt",
        minus="results/DS/{samples}/{samples}_{build}_ds_sim_minus_{regions}_combined_info.txt",
    output:
        plus="results/DS/{samples}/{samples}_{build}_ds_sim_plus_{regions}_combined_rpkm.txt",
        minus="results/DS/{samples}/{samples}_{build}_ds_sim_minus_{regions}_combined_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_{regions}_ds_sim.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_{regions}_ds_sim.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Calculating RPKM values (plus strand)..." &&
        python3 workflow/scripts/RPKM.py \
        -i {input.plus} \
        -chse 2 3 \
        -c 7 \
        -mr 0 \
        -o {output.plus} &&
        echo "`date -R`: Success! RPKMs are calculated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Calculating RPKM values (minus strand)..." &&
        python3 workflow/scripts/RPKM.py \
        -i {input.minus} \
        -chse 2 3 \
        -c 7 \
        -mr 0 \
        -o {output.minus} &&
        echo "`date -R`: Success! RPKMs are calculated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """
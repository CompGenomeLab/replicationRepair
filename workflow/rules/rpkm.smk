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

rule rpkm_sim:
    input:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{regions}_combined_info.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{regions}_combined_info.txt",
    output:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{regions}_combined_rpkm.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{regions}_combined_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_{method}_sim_{regions}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_{method}_sim_{regions}.benchmark.txt",
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

rule rpkm_xr_intergenic:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}_intergenic_combined_info.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}_intergenic_combined_info.txt",
    output:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}_intergenic_combined_rpkm.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}_intergenic_combined_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_{regions}_xr_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_{regions}_xr_intergenic.benchmark.txt",
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

rule rpkm_ds_intergenic:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}_intergenic_combined_info.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}_intergenic_combined_info.txt",
    output:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}_intergenic_combined_rpkm.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}_intergenic_combined_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_{regions}_ds_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_{regions}_ds_intergenic.benchmark.txt",
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

rule rpkm_sim_intergenic:
    input:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{regions}_intergenic_combined_info.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{regions}_intergenic_combined_info.txt",
    output:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{regions}_intergenic_combined_rpkm.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{regions}_intergenic_combined_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_{method}_sim_{regions}_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_{method}_sim_{regions}_intergenic.benchmark.txt",
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

rule rpkm_noWindows_xr:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{noWindows}_info.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{noWindows}_info.txt",
    output:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{noWindows}_rpkm.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{noWindows}_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_noWindows_{noWindows}_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_noWindows_{noWindows}_xr.benchmark.txt",
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

rule rpkm_noWindows_ds:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{noWindows}_info.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{noWindows}_info.txt",
    output:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{noWindows}_rpkm.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{noWindows}_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_noWindows_{noWindows}_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_noWindows_{noWindows}_ds.benchmark.txt",
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

rule rpkm_noWindows_sim:
    input:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{noWindows}_info.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{noWindows}_info.txt",
    output:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{noWindows}_rpkm.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{noWindows}_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_noWindows_{method}_sim_{regions}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_noWindows_{method}_sim_{regions}.benchmark.txt",
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

rule rpkm_noWindows_xr_intergenic:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{noWindows}_intergenic_info.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{noWindows}_intergenic_info.txt",
    output:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{noWindows}_intergenic_rpkm.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{noWindows}_intergenic_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_noWindows_{noWindows}_xr_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_noWindows_{noWindows}_xr_intergenic.benchmark.txt",
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

rule rpkm_noWindows_ds_intergenic:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{noWindows}_intergenic_info.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{noWindows}_intergenic_info.txt",
    output:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{noWindows}_intergenic_rpkm.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{noWindows}_intergenic_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_noWindows_{noWindows}_ds_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_noWindows_{noWindows}_ds_intergenic.benchmark.txt",
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

rule rpkm_noWindows_sim_intergenic:
    input:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{noWindows}_intergenic_info.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{noWindows}_intergenic_info.txt",
    output:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{noWindows}_intergenic_rpkm.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{noWindows}_intergenic_rpkm.txt",
    log:
        "logs/{samples}/{samples}_{build}_rpkm_noWindows_{method}_sim_{regions}_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_rpkm_noWindows_{method}_sim_{regions}_intergenic.benchmark.txt",
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

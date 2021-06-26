rule combine_files:
    input:
        lambda w: combineOutputs(config["build"], config["sample_xr"], config["sample_ds"], w.regions),
    output:
        "results/final_reports_{build}_{regions}.txt",
    log:
        "logs/{build}_combine_files_{regions}.log",
    benchmark:
        "logs/{build}_combine_files_{regions}.benchmark.txt",
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
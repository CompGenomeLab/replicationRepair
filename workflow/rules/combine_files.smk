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
        (echo "`date -R`: Combining files..." &&
        cat {input} > {output} &&
        echo "`date -R`: Success! Files are combined." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """
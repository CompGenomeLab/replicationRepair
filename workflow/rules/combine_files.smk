rule combine_files:
    input:
        lambda w: combineOutputs(config["build"], config["sample_xr"], config["sample_ds"], w.regions),
    output:
        "results/final/final_reports_{build}_{regions}.txt",
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

rule combine_files_sim:
    input:
        lambda w: combineOutputs(config["build"], config["sample_xr"], config["sample_ds"], w.regions, outformat="sim"),
    output:
        "results/final/final_reports_sim_{build}_{regions}.txt",
    log:
        "logs/{build}_combine_files_sim_{regions}.log",
    benchmark:
        "logs/{build}_combine_files_sim_{regions}.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combining files..." &&
        cat {input} > {output} &&
        echo "`date -R`: Success! Files are combined." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """

rule combine_files_intergenic:
    input:
        lambda w: combineOutputs(config["build"], config["sample_xr"], config["sample_ds"], w.regions, outformat="intergenic"),
    output:
        "results/final/final_reports_{build}_{regions}_intergenic.txt",
    log:
        "logs/{build}_combine_files_{regions}_intergenic.log",
    benchmark:
        "logs/{build}_combine_files_{regions}_intergenic.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combining files..." &&
        cat {input} > {output} &&
        echo "`date -R`: Success! Files are combined." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """

rule combine_files_sim_intergenic:
    input:
        lambda w: combineOutputs(config["build"], config["sample_xr"], config["sample_ds"], w.regions, outformat="sim_intergenic"),
    output:
        "results/final/final_reports_sim_{build}_{regions}_intergenic.txt",
    log:
        "logs/{build}_combine_files_sim_{regions}_intergenic.log",
    benchmark:
        "logs/{build}_combine_files_sim_{regions}_intergenic.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combining files..." &&
        cat {input} > {output} &&
        echo "`date -R`: Success! Files are combined." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """
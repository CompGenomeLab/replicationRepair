rule combine_files:
    input:
        lambda w: combineOutputs(config["xr"]["samples"], config["ds"]["samples"], w.regions),
    output:
        "results/final/final_reports_hg19_{regions}.txt",
    log:
        "logs/rule/analysis/hg19_combine_files_{regions}.log",
    benchmark:
        "logs/rule/analysis/hg19_combine_files_{regions}.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combining files..." &&
        cat {input} > {output} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule combine_files_sim:
    input:
        lambda w: combineOutputs(config["xr"]["samples"], config["ds"]["samples"], w.regions, outformat="sim"),
    output:
        "results/final/final_reports_sim_hg19_{regions}.txt",
    log:
        "logs/rule/analysis/hg19_combine_files_sim_{regions}.log",
    benchmark:
        "logs/rule/analysis/hg19_combine_files_sim_{regions}.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combining files..." &&
        cat {input} > {output} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule combine_files_intergenic:
    input:
        lambda w: combineOutputs(config["xr"]["samples"], config["ds"]["samples"], w.regions, outformat="intergenic"),
    output:
        "results/final/final_reports_hg19_{regions}_intergenic.txt",
    log:
        "logs/rule/analysis/hg19_combine_files_{regions}_intergenic.log",
    benchmark:
        "logs/rule/analysis/hg19_combine_files_{regions}_intergenic.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combining files..." &&
        cat {input} > {output} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

rule combine_files_sim_intergenic:
    input:
        lambda w: combineOutputs(config["xr"]["samples"], config["ds"]["samples"], w.regions, outformat="sim_intergenic"),
    output:
        "results/final/final_reports_sim_hg19_{regions}_intergenic.txt",
    log:
        "logs/rule/analysis/hg19_combine_files_sim_{regions}_intergenic.log",
    benchmark:
        "logs/rule/analysis/hg19_combine_files_sim_{regions}_intergenic.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combining files..." &&
        cat {input} > {output} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

rule combine_files_tss:
    input:
        real=lambda w: combineOutputs(config["xr"]["samples"], config["ds"]["samples"], outformat=w.tss_tes),
        sim=lambda w: combineOutputs(config["xr"]["samples"], config["ds"]["samples"], outformat=(w.tss_tes + "_sim")),
    output:
        real="results/final/final_reports_hg19_{tss_tes}.txt",
        sim="results/final/final_reports_sim_hg19_{tss_tes}.txt",
    log:
        "logs/rule/analysis/hg19_combine_files_{tss_tes}.log",
    benchmark:
        "logs/rule/analysis/hg19_combine_files_{tss_tes}.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combining files..." &&
        cat {input.real} > {output.real} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Combining simulated files..." &&
        cat {input.sim} > {output.sim} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
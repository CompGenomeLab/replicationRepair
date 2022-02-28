rule combine_windows_mutation:
    input:
        "results/mutation/{samples}/{samples}_target_mut_comb_{regions}_org.txt",
    output:
        "results/mutation/{samples}/{samples}_target_mut_comb_{regions}_combined.txt",
    params:
        lambda w: getCombine(w.regions, config["region_mut_comb_opt"], config["regions_mut"]), 
    log:
        "logs/rule/analysis/{samples}/{samples}_combine_windows_{regions}_mutation.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_combine_windows_{regions}_mutation.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combine windows..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {input} \
        {params} \
        -o {output} &&
        echo "`date -R`: Success! Windows are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

rule combine_windows_mutation_intergenic:
    input:
        "results/mutation/{samples}/{samples}_target_mut_comb_{regions}_intergenic_org.txt",
    output:
        "results/mutation/{samples}/{samples}_target_mut_comb_{regions}_intergenic_combined.txt",
    params:
        lambda w: getCombine(w.regions, config["region_mut_comb_opt"], config["regions_mut"]), 
    log:
        "logs/rule/analysis/{samples}/{samples}_combine_windows_{regions}_mutation_intergenic.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_combine_windows_{regions}_mutation_intergenic.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Combine windows..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {input} \
        {params} \
        -o {output} &&
        echo "`date -R`: Success! Windows are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """
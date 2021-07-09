rule intersect_xr:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_damSite.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_damSite.txt",
    output:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}.txt",
        region=temp("results/{samples}_{build}_{regions}_xr.bed"),
    params:
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
    log:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_xr.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """   
        (echo "`date -R`: Copy region file {params.region}..." &&
        cp {params.region} {output.region} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1   
          
        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intersect_ds:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_damSite.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_damSite.txt",
    output:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}.txt",
        region=temp("results/{samples}_{build}_{regions}_ds.bed"),
    params:
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),    
    log:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_ds.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """   
        (echo "`date -R`: Copy region file {params.region}..." &&
        cp {params.region} {output.region} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1   
          
        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intersect_sim:
    input:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_damSite.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_damSite.txt",
    output:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{regions}.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{regions}.txt",
        region=temp("results/{samples}_{build}_{method}_{regions}_sim.bed"),
    params:
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
    log:
        "logs/{samples}/{samples}_{build}_intersect_sim_{regions}_{method}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intersect_sim_{regions}_{method}.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """   
        (echo "`date -R`: Copy region file {params.region}..." &&
        cp {params.region} {output.region} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1   
          
        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intersect_xr_intergenic:
    input:
        plus="results/XR/{samples}/{samples}_{build}_intergenic_sorted_plus.bed",
        minus="results/XR/{samples}/{samples}_{build}_intergenic_sorted_minus.bed",
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
    output:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}_intergenic.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}_intergenic.txt",
        region=temp("results/{samples}_{build}_{regions}_xr_intergenic.bed"),
    log:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_xr_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_xr_intergenic.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """   
        (echo "`date -R`: Copy region file {input.region}..." &&
        cp {input.region} {output.region} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1   
          
        (echo "`date -R`: Intersecting plus strand with {input.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {input.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intersect_ds_intergenic:
    input:
        plus="results/DS/{samples}/{samples}_{build}_intergenic_sorted_ds_dipyrimidines_plus.bed",
        minus="results/DS/{samples}/{samples}_{build}_intergenic_sorted_ds_dipyrimidines_minus.bed",
    output:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}_intergenic.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}_intergenic.txt",
        region=temp("results/{samples}_{build}_{regions}_ds_intergenic.bed"),
    params:
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),    
    log:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_ds_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_ds_intergenic.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """   
        (echo "`date -R`: Copy region file {params.region}..." &&
        cp {params.region} {output.region} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1   
          
        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intersect_sim_intergenic:
    input:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_intergenic.bed",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_intergenic.bed",
    output:
        plus="results/sim/{samples}/{samples}_{build}_{method}_sim_plus_{regions}_intergenic.txt",
        minus="results/sim/{samples}/{samples}_{build}_{method}_sim_minus_{regions}_intergenic.txt",
        region=temp("results/{samples}_{build}_{method}_{regions}_sim_intergenic.bed"),
    params:
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
    log:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_{method}_sim_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_{method}_sim_intergenic.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """   
        (echo "`date -R`: Copy region file {params.region}..." &&
        cp {params.region} {output.region} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1   
          
        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intersect_mutation:
    input:
        plus="results/mutation/{samples}/{samples}_target_mut_plus.tsv",
        minus="results/mutation/{samples}/{samples}_target_mut_minus.tsv",
    output:
        plus="results/mutation/{samples}/{samples}_target_mut_plus_{regions}.txt",
        minus="results/mutation/{samples}/{samples}_target_mut_minus_{regions}.txt",
        region=temp("results/{samples}_{regions}_mut.bed"),
    params:
        region=lambda w: getRegion(w.regions, config["region_mut_file"], config["regions_mut"]),    
    log:
        "logs/{samples}/{samples}_intersect_mutation_{regions}.log",
    benchmark:
        "logs/{samples}/{samples}_intersect_mutation_{regions}.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """   
        (echo "`date -R`: Copy region file {params.region}..." &&
        cp {params.region} {output.region} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1   
          
        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intersect_mutation_intergenic:
    input:
        plus="results/mutation/{samples}/{samples}_intergenic_target_mut_plus.tsv",
        minus="results/mutation/{samples}/{samples}_intergenic_target_mut_minus.tsv",
    output:
        plus="results/mutation/{samples}/{samples}_target_mut_plus_{regions}_intergenic.txt",
        minus="results/mutation/{samples}/{samples}_target_mut_minus_{regions}_intergenic.txt",
        region=temp("results/{samples}_{regions}_mut_intergenic.bed"),
    params:
        region=lambda w: getRegion(w.regions, config["region_mut_file"], config["regions_mut"]),    
    log:
        "logs/{samples}/{samples}_intersect_mutation_{regions}_intergenic.log",
    benchmark:
        "logs/{samples}/{samples}_intersect_mutation_{regions}_intergenic.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """   
        (echo "`date -R`: Copy region file {params.region}..." &&
        cp {params.region} {output.region} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1   
          
        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """
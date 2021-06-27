rule intersect_xr:
    input:
        plus="resources/samples/XR/{samples}_{build}_sorted_plus.bed",
        minus="resources/samples/XR/{samples}_{build}_sorted_minus.bed",
    output:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}.txt",
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
        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -sorted -a results/regions/{params.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -sorted -a results/regions/{params.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intersect_ds:
    input:
        plus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
    output:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}.txt",
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
        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -sorted -a results/regions/{params.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -sorted -a results/regions/{params.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intersect_xr_sim:
    input:
        "resources/samples/XR/{samples}_{build}_xr_sim.bed",
    output:
        plus_sim=temp("results/XR/{samples}/{samples}_{build}_xr_sim_plus.bed"),
        minus_sim=temp("results/XR/{samples}/{samples}_{build}_xr_sim_minus.bed"),
        plus="results/XR/{samples}/{samples}_{build}_xr_sim_plus_{regions}.txt",
        minus="results/XR/{samples}/{samples}_{build}_xr_sim_minus_{regions}.txt",
    params:
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
    log:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_xr_sim.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_xr_sim.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Separating plus stranded reads..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus_sim} &&
        echo "`date -R`: Success! Reads are separated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Separating minus stranded reads..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus_sim} &&
        echo "`date -R`: Success! Reads are separated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -a results/regions/{params.region} \
        -b {output.plus_sim} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -a results/regions/{params.region} \
        -b {output.minus_sim} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intersect_ds_sim:
    input:
        "resources/samples/DS/{samples}_{build}_ds_sim.bed",
    output:
        plus_sim=temp("results/DS/{samples}/{samples}_{build}_ds_sim_plus.bed"),
        minus_sim=temp("results/DS/{samples}/{samples}_{build}_ds_sim_minus.bed"),
        plus="results/DS/{samples}/{samples}_{build}_ds_sim_plus_{regions}.txt",
        minus="results/DS/{samples}/{samples}_{build}_ds_sim_minus_{regions}.txt",
    params:
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),    
    log:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_ds_sim.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_ds_sim.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Separating plus stranded reads..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus_sim} &&
        echo "`date -R`: Success! Reads are separated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Separating minus stranded reads..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus_sim} &&
        echo "`date -R`: Success! Reads are separated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -a results/regions/{params.region} \
        -b {output.plus_sim} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -a results/regions/{params.region} \
        -b {output.minus_sim} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """


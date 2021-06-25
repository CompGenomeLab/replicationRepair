rule intersect_xr:
    input:
        plus="resources/samples/XR/{samples}_{build}_sorted_plus.bed",
        minus="resources/samples/XR/{samples}_{build}_sorted_minus.bed",
        region=expand("results/regions/{region}", region = config["regions"]),
    output:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_{regions}.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_{regions}.txt",
    log:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_xr.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Intersecting plus strand with {input.region}..." &&
        bedtools intersect \
        -sorted -a {input.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {input.region}..." &&
        bedtools intersect \
        -sorted -a {input.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule intersect_ds:
    input:
        plus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
        region=expand("results/regions/{region}", region = config["regions"]),
    output:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_{regions}.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_{regions}.txt",
    log:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_intersect_{regions}_ds.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Intersecting plus strand with {input.region}..." &&
        bedtools intersect \
        -sorted -a {input.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {input.region}..." &&
        bedtools intersect \
        -sorted -a {input.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """
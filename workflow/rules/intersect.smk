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
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1   
          
        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
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
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1   
          
        (echo "`date -R`: Intersecting plus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.plus} \
        -wa -c -F 0.5 > {output.plus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Intersecting minus strand with {params.region}..." &&
        bedtools intersect \
        -a {output.region} \
        -b {input.minus} \
        -wa -c -F 0.5 > {output.minus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
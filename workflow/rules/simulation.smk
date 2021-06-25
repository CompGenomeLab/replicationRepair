
rule simulation:
    input:
        bed="results/input/{samples}/{samples}_{build}_sorted_chr.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
        regions="resources/ref_genomes/{build}/genome_{build}.ron", 
    output:
        fa=temp("results/input/{samples}/{samples}_{build}_sorted_chr.fa"),
        simbed="results/input/{samples}/{samples}_{build}_sim.bed",
        sim="results/input/{samples}/{samples}_{build}_sim.fa",
    log:
        "logs/{samples}/{samples}_{build}_simulation.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_simulation.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output.fa} \
        -s &&
        echo "`date -R`: Success! {input.bed} is converted." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Simulating reads..." &&
        boquila \
        --fasta {output.fa} \
        --bed {output.simbed} \
        --ref {input.genome} \
        --seed 1 \
        --regions {input.regions} \
        > {output.sim} &&
        echo "`date -R`: Success! Simulation is done." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """
rule simulation_ds_input:
    input:
        plus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus_{kmer}.bed",
        minus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus_{kmer}.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
        inpfile=lambda w: getInput(w.samples, config["input"]["exist"], config["input"]["files"], config["input"]["sample"], config["ds"]["samples"], config["build"], region="IZ"),
    output:
        bed=temp("resources/samples/sim/{samples}_{build}_sorted_ds_{kmer}_dipyrimidines.bed"),
        fa=temp("resources/samples/sim/{samples}_{build}_sorted_ds_{kmer}_dipyrimidines.fa"),
        sim="resources/samples/sim/{samples}_{build}_ds_sim_{kmer}.fa",
        sam=temp("resources/samples/sim/{samples}_cutadapt_ds_sim_{kmer}_{build}.sam"),
        bam=temp("resources/samples/sim/{samples}_cutadapt_ds_sim_{kmer}_{build}.bam"),
        bam_sorted=temp("resources/samples/sim/{samples}_cutadapt_ds_sim_{kmer}_{build}_sorted.bam"),
        idx=temp("resources/samples/sim/{samples}_cutadapt_ds_sim_{kmer}_{build}_sorted.bam.bai"),
        simbed="resources/samples/sim/{samples}_{build}_ds_sim_{kmer}.bed",  
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        kmer="1",
    log:
        "logs/{samples}/{samples}_{build}_simulation_{kmer}_ds_input.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_simulation_{kmer}_ds_input.benchmark.txt",
    conda:
        "../envs/simulation.yaml"
    shell:
        """
        (echo "`date -R`: Combine files..." &&
        cat {input.plus} {input.minus} > {output.bed} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        
        (echo "`date -R`: Converting {output.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {output.bed} \
        -fo {output.fa} \
        -s &&
        echo "`date -R`: Success! {output.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Simulating reads with input file..." &&
        /cta/users/uakkose/bin/boquila_cem_skip \
        --fasta {output.fa} \
        --inseqFasta \
        --inseq {input.inpfile} \
        --kmer {params.kmer} \
        --seed 1 \
        > {output.sim} &&
        echo "`date -R`: Success! Simulation is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Aligning fasta file..." &&
        bowtie2 \
        -x {params.ref_genome} \
        -f {output.sim} -S {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools view -Sbh -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Sorting (coordinates) bam file..." &&
        samtools sort {output.bam} > {output.bam_sorted} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Index bam file..." &&
        samtools index {output.bam_sorted} {output.idx} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." && 
        samtools view -b {output.bam_sorted} |&
        bedtools bamtobed > {output.simbed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output.simbed}; exit 1; }}  ) >> {log} 2>&1
        """

rule simulation_xr_input:
    input:
        plus="resources/samples/XR/{samples}_{build}_sorted_plus_{kmer}.bed",
        minus="resources/samples/XR/{samples}_{build}_sorted_minus_{kmer}.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2", 
        inpfile=lambda w: getInput(w.samples, config["input"]["exist"], config["input"]["files"], config["input"]["sample"], config["xr"]["samples"], config["build"], region="IZ"),
    output:
        bed=temp("resources/samples/sim/{samples}_{build}_sorted_chr_{kmer}.bed"),
        fa=temp("resources/samples/sim/{samples}_{build}_sorted_chr_{kmer}.fa"),
        sim="resources/samples/sim/{samples}_{build}_xr_sim_{kmer}.fa",
        sam=temp("resources/samples/sim/{samples}_cutadapt_xr_sim_{kmer}_{build}.sam"),
        bam=temp("resources/samples/sim/{samples}_cutadapt_xr_sim_{kmer}_{build}.bam"),
        bam_sorted=temp("resources/samples/sim/{samples}_cutadapt_xr_sim_{kmer}_{build}_sorted.bam"),
        idx=temp("resources/samples/sim/{samples}_cutadapt_xr_sim_{kmer}_{build}_sorted.bam.bai"),
        simbed="resources/samples/sim/{samples}_{build}_xr_sim_{kmer}.bed",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        kmer="1",
    log:
        "logs/{samples}/{samples}_{build}_simulation_{kmer}_xr_input.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_simulation_{kmer}_xr_input.benchmark.txt",
    conda:
        "../envs/simulation.yaml"
    shell:
        """
        (echo "`date -R`: Combine files..." &&
        cat {input.plus} {input.minus} > {output.bed} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Converting {output.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {output.bed} \
        -fo {output.fa} \
        -s &&
        echo "`date -R`: Success! {output.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Simulating reads with input file..." &&
        /cta/users/uakkose/bin/boquila_cem_skip \
        --fasta {output.fa} \
        --inseqFasta \
        --inseq {input.inpfile} \
        --kmer {params.kmer} \
        --seed 1 \
        > {output.sim} &&
        echo "`date -R`: Success! Simulation is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Aligning fasta file..." &&
        bowtie2 \
        -x {params.ref_genome} \
        -f {output.sim} -S {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools view -Sbh -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Sorting (coordinates) bam file..." &&
        samtools sort {output.bam} > {output.bam_sorted} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Index bam file..." &&
        samtools index {output.bam_sorted} {output.idx} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." && 
        samtools view -b {output.bam_sorted} |&
        bedtools bamtobed > {output.simbed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output.simbed}; exit 1; }}  ) >> {log} 2>&1
        """


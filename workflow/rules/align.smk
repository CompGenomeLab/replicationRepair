rule bowtie2_se_okseq:
    input:
        sample=["results/okseq/{samples}/{samples}_cutadapt.fastq.gz"],
        bowtie2="resources/ref_genomes/hg19/Bowtie2/genome_hg19.1.bt2",
    output:
        sam=temp("results/okseq/{samples}/{samples}_se_hg19.sam"),
        bam=temp("results/okseq/{samples}/{samples}_se_hg19.bam"),
        bam_filt=temp("results/okseq/{samples}/{samples}_se_hg19_filt.bam"),
        bam_org="results/okseq/{samples}/{samples}_se_hg19_sorted.bam",
        idx="results/okseq/{samples}/{samples}_se_hg19_sorted.bam.bai",
    params:
        ref_genome="resources/ref_genomes/hg19/Bowtie2/genome_hg19",
        extra="--seed 1",
    threads: 16  
    log:
        "logs/rule/analysis/{samples}/{samples}_hg19_bowtie2_se_okseq.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_hg19_bowtie2_se_okseq.benchmark.txt",
    conda:
        "../envs/align.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq file..." &&
        bowtie2 \
        --threads {threads} \
        {params.extra} \
        -x {params.ref_genome} \
        -U {input.sample[0]} -S {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools view -q 20 -Sbh -o {output.bam} {output.sam} &&
        samtools rmdup -s {output.bam} {output.bam_filt} &&
        samtools sort {output.bam_filt} -o {output.bam_org}
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        
        (echo "`date -R`: Indexing bam..." &&
        samtools index {output.bam_org}
        echo "`date -R`: Success! Indexing is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """

rule bowtie2_se_edu:
    input:
        sample=["resources/samples/edu/{samples}.fastq.gz"],
        bowtie2="resources/ref_genomes/hg19/Bowtie2/genome_hg19.1.bt2",
    output:
        sam=temp("results/edu/{samples}/{samples}_se_hg19.sam"),
        bam="results/edu/{samples}/{samples}_se_hg19.bam",
    params:
        ref_genome="resources/ref_genomes/hg19/Bowtie2/genome_hg19",
        extra="--seed 1",
    threads: 16  
    log:
        "logs/rule/analysis/{samples}/{samples}_hg19_bowtie2_se_edu.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_hg19_bowtie2_se_edu.benchmark.txt",
    conda:
        "../envs/align.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq file..." &&
        bowtie2 \
        --threads {threads} \
        {params.extra} \
        -x {params.ref_genome} \
        --no-mixed --no-discordant --reorder \
        -U {input.sample[0]} -S {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools view -Sbh -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """

rule bowtie2_pe_edu:
    input:
        sample=["resources/samples/edu/{samples}_R1.fastq.gz", "resources/samples/edu/{samples}_R2.fastq.gz"],
        bowtie2="resources/ref_genomes/hg19/Bowtie2/genome_hg19.1.bt2",
    output:
        sam=temp("results/edu/{samples}/{samples}_pe_hg19.sam"),
        bam="results/edu/{samples}/{samples}_pe_hg19.bam",
    params:
        ref_genome="resources/ref_genomes/hg19/Bowtie2/genome_hg19",
        extra="-X 1000 --seed 1",
    threads: 16  
    log:
        "logs/rule/analysis/{samples}/{samples}_hg19_bowtie2_pe_edu.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_hg19_bowtie2_pe_edu.benchmark.txt",
    conda:
        "../envs/align.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq files..." &&
        bowtie2 \
        --threads {threads} \
        {params.extra} \
        -x {params.ref_genome} \
        --no-mixed --no-discordant --reorder -X 1000 \
        -1 {input.sample[0]} -2 {input.sample[1]} -S {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools view -Sb -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """

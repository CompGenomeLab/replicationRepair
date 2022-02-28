rule bam2bed_se_edu:
    input:
        "results/edu/{samples}/{samples}_se_{build}.bam",
    output:
        bam=temp("results/edu/{samples}/{samples}_{build}_collated.bam"),
        fix=temp("results/edu/{samples}/{samples}_{build}_collated_fixmate.bam"),
        sort=temp("results/edu/{samples}/{samples}_{build}_sorted.bam"),
        rmdup=temp("results/edu/{samples}/{samples}_{build}_sorted_rmdup.bam"),
        bed="results/edu/{samples}/{samples}_{build}_se.bed",
    params:
        q_trim="-q 20", 
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_bam2bed_se_edu.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_bam2bed_se_edu.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Collating bam file..." &&
        samtools collate -o {output.bam} {input}  &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Fixmate..." &&
        samtools fixmate -m {output.bam} {output.fix} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Sorting..." &&
        samtools sort -o {output.sort} {output.fix} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Removing duplicates..." &&
        samtools markdup -r {output.sort} {output.rmdup} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." && 
        samtools view {params.q_trim} {output.rmdup} |&
        bedtools bamtobed > {output.bed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule bam2bed_pe_edu:
    input:
        "results/edu/{samples}/{samples}_pe_{build}.bam",
    output:
        bam=temp("results/edu/{samples}/{samples}_{build}_collated.bam"),
        fix=temp("results/edu/{samples}/{samples}_{build}_collated_fixmate.bam"),
        sort=temp("results/edu/{samples}/{samples}_{build}_sorted.bam"),
        rmdup="results/edu/{samples}/{samples}_{build}_sorted_rmdup.bam",
        bed="results/edu/{samples}/{samples}_{build}_pe.bed",
    params:
        q_trim="-q 20 -bf 0x2",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_bam2bed_pe_edu.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_bam2bed_pe_edu.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell: 
        """
        (echo "`date -R`: Collating bam file..." &&
        samtools collate -o {output.bam} {input}  &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Fixmate..." &&
        samtools fixmate -m {output.bam} {output.fix} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Sorting..." &&
        samtools sort -o {output.sort} {output.fix} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Removing duplicates..." &&
        samtools markdup -r {output.sort} {output.rmdup} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Indexing..." &&
        samtools index {output.rmdup} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." &&
        samtools view {params.q_trim} {output.rmdup} |&
        bedtools bamtobed -bedpe -mate1 |&
        awk '{{\
            if ($9=="+")\
                print $1"\\t"$2"\\t"$6"\\t"$7"\\t"$8"\\t"$9;\
            else if ($9=="-")\
                print $1"\\t"$5"\\t"$3"\\t"$7"\\t"$8"\\t"$9;\
            }}' > {output.bed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """

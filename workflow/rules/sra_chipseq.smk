rule sra_se_chipseq:
    output:
        "resources/samples/chipseq/{samples}.fastq.gz", 
    params:
        srr=lambda w: getSRR(w.samples, config["chipseq"]["srr"]["codes"], 
            config["chipseq"]["samples"]),
        name="{samples}",
    log:
        "logs/{samples}/{samples}_se_sra_chipseq.log",
    benchmark:
        "logs/{samples}/{samples}_se_sra_chipseq.benchmark.txt",
    wildcard_constraints:
        samples='|'.join([x for x in config["chipseq"]["samples"]])
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/samples/chipseq/{params.name}.fastq
        touch {log}
        echo "`date -R`: single-end layout" >> {log}

        srrList=$(echo {params.srr} | tr ":" "\\n")
        echo $srrList >> {log}

        for srr in $srrList; do
            (echo "`date -R`: Downloading $srr files..." &&
            fasterq-dump \
            --threads {threads} \
            --progress $srr \
            -t resources/samples/chipseq/ \
            -o resources/samples/chipseq/${{srr}}.fastq &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
            >> {log} 2>&1 

            cat resources/samples/chipseq/${{srr}}.fastq >> resources/samples/chipseq/{params.name}.fastq

            rm resources/samples/chipseq/${{srr}}.fastq ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/samples/chipseq/{params.name}.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
        >> {log} 2>&1
        """

rule sra_pe_chipseq:
    output:
        "resources/samples/chipseq/{samples}_1.fastq.gz", 
        "resources/samples/chipseq/{samples}_2.fastq.gz", 
    params:
        srr=lambda w: getSRR(w.samples, config["chipseq"]["srr"]["codes"], 
            config["chipseq"]["samples"]),
        name="{samples}",
    log:
        "logs/{samples}/{samples}_pe_sra_chipseq.log",
    benchmark:
        "logs/{samples}/{samples}_pe_sra_chipseq.benchmark.txt",
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/samples/chipseq/{params.name}_1.fastq
        touch resources/samples/chipseq/{params.name}_2.fastq
        touch {log}
        echo "`date -R`: paired-end layout" >> {log}

        srrList=$(echo {params.srr} | tr ":" "\\n")
        echo $srrList >> {log}

        for srr in $srrList; do
            (echo "`date -R`: Downloading $srr files..." &&
            fasterq-dump \
            --threads {threads} \
            --progress $srr \
            --split-files \
            -O resources/samples/chipseq/ &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
            >> {log} 2>&1

            (cat resources/samples/chipseq/${{srr}}_1.fastq >> resources/samples/chipseq/{params.name}_1.fastq) >> {log} 2>&1
            (cat resources/samples/chipseq/${{srr}}_2.fastq >> resources/samples/chipseq/{params.name}_2.fastq) >> {log} 2>&1

            (rm resources/samples/chipseq/${{srr}}_1.fastq) >> {log} 2>&1
            (rm resources/samples/chipseq/${{srr}}_2.fastq) >> {log} 2>&1 ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/samples/chipseq/{params.name}_1.fastq &&
        gzip resources/samples/chipseq/{params.name}_2.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
        >> {log} 2>&1
        """

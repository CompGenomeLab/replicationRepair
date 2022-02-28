rule sra_se_okseq:
    output:
        "resources/samples/okseq/{samples}.fastq.gz", 
    params:
        srr=lambda w: getSRR(w.samples, config["okseq"]["srr"]["codes"], 
            config["okseq"]["samples"]),
        name="{samples}",
    log:
        "logs/rule/analysis/{samples}/{samples}_se_sra_okseq.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_se_sra_okseq.benchmark.txt",
    wildcard_constraints:
        samples='|'.join([x for x in config["okseq"]["samples"]])
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/samples/okseq/{params.name}.fastq
        touch {log}
        echo "`date -R`: single-end layout" >> {log}

        srrList=$(echo {params.srr} | tr ":" "\\n")
        echo $srrList >> {log}

        for srr in $srrList; do
            (echo "`date -R`: Downloading $srr files..." &&
            fasterq-dump \
            --threads {threads} \
            --progress $srr \
            -t resources/samples/okseq/ \
            -o resources/samples/okseq/${{srr}}.fastq &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
            >> {log} 2>&1 

            cat resources/samples/okseq/${{srr}}.fastq >> resources/samples/okseq/{params.name}.fastq

            rm resources/samples/okseq/${{srr}}.fastq ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/samples/okseq/{params.name}.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
        >> {log} 2>&1
        """

rule sra_pe_okseq:
    output:
        "resources/samples/okseq/{samples}_1.fastq.gz", 
        "resources/samples/okseq/{samples}_2.fastq.gz", 
    params:
        srr=lambda w: getSRR(w.samples, config["okseq"]["srr"]["codes"], 
            config["okseq"]["samples"]),
        name="{samples}",
    log:
        "logs/rule/analysis/{samples}/{samples}_pe_sra_okseq.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_pe_sra_okseq.benchmark.txt",
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/samples/okseq/{params.name}_1.fastq
        touch resources/samples/okseq/{params.name}_2.fastq
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
            -O resources/samples/okseq/ &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
            >> {log} 2>&1

            (cat resources/samples/okseq/${{srr}}_1.fastq >> resources/samples/okseq/{params.name}_1.fastq) >> {log} 2>&1
            (cat resources/samples/okseq/${{srr}}_2.fastq >> resources/samples/okseq/{params.name}_2.fastq) >> {log} 2>&1

            (rm resources/samples/okseq/${{srr}}_1.fastq) >> {log} 2>&1
            (rm resources/samples/okseq/${{srr}}_2.fastq) >> {log} 2>&1 ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/samples/okseq/{params.name}_1.fastq &&
        gzip resources/samples/okseq/{params.name}_2.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
        >> {log} 2>&1
        """

rule sra_se_edu:
    output:
        "resources/samples/edu/{samples}.fastq.gz", 
    params:
        srr=lambda w: getSRR(w.samples, config["edu"]["srr"]["codes"], 
            config["edu"]["samples"]),
        name="{samples}",
    log:
        "logs/{samples}/{samples}_se_sra_edu.log",
    benchmark:
        "logs/{samples}/{samples}_se_sra_edu.benchmark.txt",
    wildcard_constraints:
        samples='|'.join([x for x in config["edu"]["samples"]])
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/samples/edu/{params.name}.fastq
        touch {log}
        echo "`date -R`: single-end layout" >> {log}

        srrList=$(echo {params.srr} | tr ":" "\\n")
        echo $srrList >> {log}

        for srr in $srrList; do
            (echo "`date -R`: Downloading $srr files..." &&
            fasterq-dump \
            --threads {threads} \
            --progress $srr \
            -t resources/samples/edu/ \
            -o resources/samples/edu/${{srr}}.fastq &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
            >> {log} 2>&1 

            cat resources/samples/edu/${{srr}}.fastq >> resources/samples/edu/{params.name}.fastq

            rm resources/samples/edu/${{srr}}.fastq ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/samples/edu/{params.name}.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
        >> {log} 2>&1
        """

rule sra_pe_edu:
    output:
        "resources/samples/edu/{samples}_1.fastq.gz", 
        "resources/samples/edu/{samples}_2.fastq.gz", 
    params:
        srr=lambda w: getSRR(w.samples, config["edu"]["srr"]["codes"], 
            config["edu"]["samples"]),
        name="{samples}",
    log:
        "logs/{samples}/{samples}_pe_sra_edu.log",
    benchmark:
        "logs/{samples}/{samples}_pe_sra_edu.benchmark.txt",
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/samples/edu/{params.name}_1.fastq
        touch resources/samples/edu/{params.name}_2.fastq
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
            -O resources/samples/edu/ &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
            >> {log} 2>&1

            (cat resources/samples/edu/${{srr}}_1.fastq >> resources/samples/edu/{params.name}_1.fastq) >> {log} 2>&1
            (cat resources/samples/edu/${{srr}}_2.fastq >> resources/samples/edu/{params.name}_2.fastq) >> {log} 2>&1

            (rm resources/samples/edu/${{srr}}_1.fastq) >> {log} 2>&1
            (rm resources/samples/edu/${{srr}}_2.fastq) >> {log} 2>&1 ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/samples/edu/{params.name}_1.fastq &&
        gzip resources/samples/edu/{params.name}_2.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
        >> {log} 2>&1
        """

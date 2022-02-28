rule getPublicRepli:
    output:
        g1="resources/samples/edu/repliseqG1b.bam",
        s1="resources/samples/edu/repliseqS1.bam",
        s2="resources/samples/edu/repliseqS2.bam",
        s3="resources/samples/edu/repliseqS3.bam",
        s4="resources/samples/edu/repliseqS4.bam",
        g2="resources/samples/edu/repliseqG2.bam",
    log:
        "logs/rule/analysis/getPublicRepli.log",
    benchmark:
        "logs/rule/analysis/getPublicRepli.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Download {output.g1}..." &&
        wget https://www.encodeproject.org/files/ENCFF001GOF/@@download/ENCFF001GOF.bam -O {output.g1} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Index {output.g1}..." &&
        samtools index {output.g1} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Download {output.s1}..." &&
        wget https://www.encodeproject.org/files/ENCFF001GOL/@@download/ENCFF001GOL.bam -O {output.s1} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Index {output.s1}..." &&
        samtools index {output.s1} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Download {output.s2}..." &&
        wget https://www.encodeproject.org/files/ENCFF001GOP/@@download/ENCFF001GOP.bam -O {output.s2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Index {output.s2}..." &&
        samtools index {output.s2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Download {output.s3}..." &&
        wget https://www.encodeproject.org/files/ENCFF001GOR/@@download/ENCFF001GOR.bam -O {output.s3} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Index {output.s3}..." &&
        samtools index {output.s3} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Download {output.s4}..." &&
        wget https://www.encodeproject.org/files/ENCFF001GOY/@@download/ENCFF001GOY.bam -O {output.s4} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Index {output.s4}..." &&
        samtools index {output.s4} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Download {output.g2}..." &&
        wget https://www.encodeproject.org/files/ENCFF001GOI/@@download/ENCFF001GOI.bam -O {output.g2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: Index {output.g2}..." &&
        samtools index {output.g2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """

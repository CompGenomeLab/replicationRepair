rule bdg2bed:
    input:
        "results/edu/early_late.bdg",    
    output:
        "results/regions/early_late.bed",
    log:
        "logs/bdg2bed/bdg2bed.log",
    benchmark:
        "logs/bdg2bed/bdg2bed.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Converting bedgraph to bed..." &&
        awk '{{print $1"\\t"$2"\\t"$3"\\t"".""\\t"$4"\\t""."}}' {input} > {output} &&
        echo "`date -R`: Success! Bed file is created." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

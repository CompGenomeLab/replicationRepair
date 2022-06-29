rule merge_bdg:
    input:
        "results/edu/merge_uv_RT.txt",
    output:
        "results/regions/merge_uv_RT.bed",
    log:
        "logs/merge_bdg/merge_bdg.log",
    benchmark:
        "logs/merge_bdg/merge_bdg.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Get the average of 2 replicates..." &&
        awk '{{print $1"\\t"$2"\\t"$3"\\t"".""\\t"$4+$5/2"\\t""."}}' \
        {input} \
        > {output} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """
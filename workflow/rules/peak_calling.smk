rule peak_calling_edu:
    input:
        lambda w: input4peakCalling(w, config["build"], "edu"),
    output:
        narrow="results/edu/{samples}_{build}_peaks.narrowPeak",
        broad="results/edu/{samples}_{build}_peaks.broadPeak",
    params:
        name="{samples}_{build}", 
        outdir="results/edu",
    log:
        "logs/{samples}/{samples}_{build}_peak_calling_edu.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_peak_calling_edu.benchmark.txt",
    conda:
        "../envs/peak_calling.yaml"
    shell:  
        """
        (echo "`date -R`: Peak calling (narrowPeak)..." && 
        macs2 callpeak \
        -t {input[0]} \
        -c {input[1]} \
        -f BAMPE \
        -n {params.name} \
        --outdir {params.outdir} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Peak calling (broadPeak)..." && 
        macs2 callpeak \
        -t {input[0]} \
        -c {input[1]} \
        -f BAMPE \
        -n {params.name} \
        --outdir {params.outdir} \
        --broad &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

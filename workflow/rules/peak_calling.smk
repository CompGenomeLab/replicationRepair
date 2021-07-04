rule peak_calling_edu:
    input:
        samp="results/edu/{name1}-{name2}/{name1}-{name2}_sorted.bam",
        inp="results/edu/{name1}-inp{name2}/{name1}-inp{name2}_sorted.bam",
    output:
        narrow="results/edu/{name1}-{name2}_{build}_peaks.narrowPeak",
        broad="results/edu/{name1}-{name2}_{build}_peaks.broadPeak",
    params:
        name="{name1}-{name2}_{build}_peaks", 
        outdir="results/edu",
    log:
        "logs/{name1}-{name2}/{name1}-{name2}_{build}_peak_calling_edu.log",
    benchmark:
        "logs/{name1}-{name2}/{name1}-{name2}_{build}_peak_calling_edu.benchmark.txt",
    conda:
        "../envs/peak_calling.yaml"
    shell:  
        """
        (echo "`date -R`: Peak calling (narrowPeak)..." && 
        macs2 callpeak \
        -t {input.samp} \
        -c {input.inp} \
        -f BAMPE \
        -n {params.name} \
        --outdir {params.outdir} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Peak calling (broadPeak)..." && 
        macs2 callpeak \
        -t {input.samp} \
        -c {input.inp} \
        -f BAMPE \
        -n {params.name} \
        --outdir {params.outdir} \
        --broad &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

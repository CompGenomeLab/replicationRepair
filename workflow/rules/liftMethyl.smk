rule liftMethyl:
    input:
        bed="resources/samples/methyl/{samples}_GRCh38.bed",
        ref="resources/ref_genomes/hg38ToHg19.over.chain",
    output:
        bed="results/methyl/{samples}/{samples}_GRCh38.bed",
        hg19="results/methyl/{samples}/{samples}_hg19.bed",
        unmapped="results/methyl/{samples}/{samples}_unmapped2hg19.bed",
        hg19chr="results/methyl/{samples}/{samples}_hg19_chr.bed",
    params:
        "results/methyl/{samples}"
    log:
        "logs/{samples}/{samples}_liftMethyl.log",
    benchmark:
        "logs/{samples}/{samples}_liftMethyl.benchmark.txt",
    conda:
        "../envs/liftMethyl.yaml"
    shell:  
        """
        ( echo "`date -R`: Reorganize.." &&
        cut -f1-6 {input.bed} > {output.bed} &&
        echo "`date -R`: done!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        # UCSC LiftOver
        ( echo "`date -R`: LiftOver.." &&
        workflow/scripts/liftOver \
        {output.bed} \
        {input.ref} \
        {output.hg19} \
        {output.unmapped} &&
        echo "`date -R`: done!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        # get only chromosomes
        ( echo "`date -R`: get only chromosomes.." &&
        grep -P "chr.?.?\\t" {output.hg19} > {output.hg19chr} &&
        echo "`date -R`: done!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """ 
 

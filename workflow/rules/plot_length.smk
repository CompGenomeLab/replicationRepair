
rule plot_length_input:
    input:
        "results/input/{samples}/{samples}_{build}_length_distribution.txt",
    output:
        report("results/input/{samples}/{samples}_{build}_length_distribution.png", 
                category="Length Distribution"),
    params:
        "{samples}_{build}",
    log:
        "logs/{samples}/{samples}_{build}_plot_length.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_plot_length.benchmark.txt",
    conda:
        "../envs/plot_nuc.yaml"
    shell:  
        """
        Rscript workflow/scripts/lengthDistPlot.r \
        -i {input} \
        -s {params} \
        -o {output} \
        -l {log}
        """
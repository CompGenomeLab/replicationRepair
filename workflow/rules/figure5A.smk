rule figure5A:
    input:  
        hela="results/okseq/HeLa/HeLa_hg19_HMMsegments_IZ.bed",
        gm06990="results/okseq/GM06990/GM06990_hg19_HMMsegments_IZ.bed",
        imr90="results/okseq/IMR90/IMR90_hg19_HMMsegments_IZ.bed",          
    output:
        report("results/plots/figure5A.pdf", caption="../report/figure5A.rst", category="Figures"),
    log:
        "logs/rule/fig/figure5A.log",
    benchmark:
        "logs/rule/fig/figure5A.benchmark.txt",
    conda:
        "../envs/venn.yaml",
    shell:
        """
        ( intervene venn -i {input} --names=HeLa-S3,GM06990,IMR90 -o workflow/results/plots/ ) > {log} 2>&1

        mv workflow/results/plots/Intervene_venn.pdf {output}
        """
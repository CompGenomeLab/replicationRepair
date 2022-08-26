rule figure4G:
    input:           
        df="results/final/final_reports_hg19_genome_hg19_50kb.txt",            
        counts="resources/ref_genomes/hg19/genome_hg19_50kb_counts.txt",                
    output:
        plot=report("results/plots/figure4G.pdf", caption="../report/figure4G.rst", category="Figures"),
        dfs="results/plot_dataframe/figure4_G.csv",
    log:
        "logs/rule/fig/figure4G.log",
    benchmark:
        "logs/rule/fig/figure4G.benchmark.txt",
    conda:
        "../envs/figures.yaml",
    shell:
        """
        Rscript workflow/scripts/figure4G.R \
        --df {input.df} \
        --counts {input.counts} \
        --data_prefix "results/plot_dataframe/figure4_G" \
        -o {output.plot} \
        --log {log} 
        """


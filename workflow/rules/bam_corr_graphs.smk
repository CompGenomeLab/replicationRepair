rule bam_corr_graphs:
    input:
        npz="results/edu/readCounts.npz",
    output:
        scatter="results/plots/scatterplot_PearsonCorr_bigwigScores.png",
        tab="results/edu/PearsonCorr_bigwigScores.tab",
        heatmap="results/plots/heatmap_SpearmanCorr_readCounts.png",
        tab2="results/edu/SpearmanCorr_readCounts.tab",
        pca="results/plots/PCA_readCounts.png",
    log:
        "logs/bam_corr_graphs.log",
    benchmark:
        "logs/bam_corr_graphs.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: Plotting correlation (scatter)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod pearson --skipZeros \
        --plotTitle "Pearson Correlation of Bam Files" \
        --whatToPlot scatterplot \
        -o {output.scatter} \
        --outFileCorMatrix {output.tab} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Plotting correlation (heatmap)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o {output.heatmap} \
        --outFileCorMatrix {output.tab2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: PCA analysis..." &&
        plotPCA -in {input.npz} \
        -o {output.pca} \
        -T "PCA of read counts" &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """

rule bam_corr_graphs_early:
    input:
        npz="results/edu/readCounts_early.npz",
    output:
        scatter="results/plots/scatterplot_PearsonCorr_bigwigScores_early.png",
        tab="results/edu/PearsonCorr_bigwigScores_early.tab",
        heatmap="results/plots/heatmap_SpearmanCorr_readCounts_early.png",
        tab2="results/edu/SpearmanCorr_readCounts_early.tab",
        pca="results/plots/PCA_readCounts_early.png",
    log:
        "logs/bam_corr_graphs_early.log",
    benchmark:
        "logs/bam_corr_graphs_early.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: Plotting correlation (scatter)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod pearson --skipZeros \
        --plotTitle "Pearson Correlation of Bam Files" \
        --whatToPlot scatterplot \
        -o {output.scatter} \
        --outFileCorMatrix {output.tab} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Plotting correlation (heatmap)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o {output.heatmap} \
        --outFileCorMatrix {output.tab2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: PCA analysis..." &&
        plotPCA -in {input.npz} \
        -o {output.pca} \
        -T "PCA of read counts" &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """

rule bam_corr_graphs_late:
    input:
        npz="results/edu/readCounts_late.npz",
    output:
        scatter="results/plots/scatterplot_PearsonCorr_bigwigScores_late.png",
        tab="results/edu/PearsonCorr_bigwigScores_late.tab",
        heatmap="results/plots/heatmap_SpearmanCorr_readCounts_late.png",
        tab2="results/edu/SpearmanCorr_readCounts_late.tab",
        pca="results/plots/PCA_readCounts_late.png",
    log:
        "logs/bam_corr_graphs_late.log",
    benchmark:
        "logs/bam_corr_graphs_late.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: Plotting correlation (scatter)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod pearson --skipZeros \
        --plotTitle "Pearson Correlation of Bam Files" \
        --whatToPlot scatterplot \
        -o {output.scatter} \
        --outFileCorMatrix {output.tab} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Plotting correlation (heatmap)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o {output.heatmap} \
        --outFileCorMatrix {output.tab2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: PCA analysis..." &&
        plotPCA -in {input.npz} \
        -o {output.pca} \
        -T "PCA of read counts" &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """

rule bam_corr_graphs_okseq:
    input:
        npz="results/okseq/readCounts_okseq.npz",
    output:
        scatter="results/plots/scatterplot_PearsonCorr_bigwigScores_okseq.png",
        tab="results/okseq/PearsonCorr_bigwigScores_okseq.tab",
        heatmap="results/plots/heatmap_SpearmanCorr_readCounts_okseq.png",
        tab2="results/okseq/SpearmanCorr_readCounts_okseq.tab",
        pca="results/plots/PCA_readCounts_okseq.png",
    log:
        "logs/bam_corr_graphs_okseq.log",
    benchmark:
        "logs/bam_corr_graphs_okseq.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: Plotting correlation (scatter)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod pearson --skipZeros \
        --plotTitle "Pearson Correlation of Bam Files" \
        --whatToPlot scatterplot \
        -o {output.scatter} \
        --outFileCorMatrix {output.tab} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Plotting correlation (heatmap)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o {output.heatmap} \
        --outFileCorMatrix {output.tab2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: PCA analysis..." &&
        plotPCA -in {input.npz} \
        -o {output.pca} \
        -T "PCA of read counts" &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """

rule bam_corr_graphs_chipseq:
    input:
        npz="results/chipseq/readCounts_chipseq.npz",
    output:
        scatter="results/plots/scatterplot_PearsonCorr_bigwigScores_chipseq.png",
        tab="results/chipseq/PearsonCorr_bigwigScores_chipseq.tab",
        heatmap="results/plots/heatmap_SpearmanCorr_readCounts_chipseq.png",
        tab2="results/chipseq/SpearmanCorr_readCounts_chipseq.tab",
        pca="results/plots/PCA_readCounts_chipseq.png",
    log:
        "logs/bam_corr_graphs_chipseq.log",
    benchmark:
        "logs/bam_corr_graphs_chipseq.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: Plotting correlation (scatter)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod pearson --skipZeros \
        --plotTitle "Pearson Correlation of Bam Files" \
        --whatToPlot scatterplot \
        -o {output.scatter} \
        --outFileCorMatrix {output.tab} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Plotting correlation (heatmap)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o {output.heatmap} \
        --outFileCorMatrix {output.tab2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: PCA analysis..." &&
        plotPCA -in {input.npz} \
        -o {output.pca} \
        -T "PCA of read counts" &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """

rule make_windows:
    input:
        region="results/regions/R21071354-EdUrep2-UV3-5hrls2_combined_hg19_peaks_repdomains_org.broadPeak",
        genome="resources/ref_genomes/hg19/genome_hg19.bed",
    output:
        final="results/regions/R21071354-EdUrep2-UV3-5hrls2_combined_broadpeaks_repdomains_windows_51_100.bed",
    params:
        intervalLen="100",
        windowNum="51",
        rev="",
    log:
        "logs/make_windows.log",
    benchmark:
        "logs/make_windows.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:  
        """
        (echo "`date -R`: Creating windows..." &&
        workflow/scripts/makewindows.sh \
        {input.region} \
        {input.genome} \
        {output} \
        {params.intervalLen} {params.windowNum} \
        {params.rev} &&
        echo "`date -R`: Success! Windowed file is ready." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """
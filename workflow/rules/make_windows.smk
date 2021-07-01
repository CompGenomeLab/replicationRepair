
rule make_windows:
    input:
        region="results/regions/R21061297-EdUrep-UV3-5hrls_combined_peaks_repdomains_org.broadPeak",
        genome="resources/ref_genomes/hg19/genome_hg19.bed",
    output:
        final="results/regions/R21061297-EdUrep-UV3-5hrls_combined_broadpeaks_windows_51_100.bed",
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
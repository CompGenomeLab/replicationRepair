rule getChromHMM:
    output:
        "resources/samples/chromHMM/wgEncodeAwgSegmentationChromhmmHelas3.bed",
    log:
        "logs/rule/analysis/getChromHMM.log",
    benchmark:
        "logs/rule/analysis/getChromHMM.benchmark.txt",
    shell:  
        """
        mkdir -p resources/samples/chromHMM/

        (echo "`date -R`: Downloading..." &&
        wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationChromhmmHelas3.bed.gz \
        -O {output}.gz &&
        gunzip {output}.gz &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

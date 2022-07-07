rule tss_xr:
    input:
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes_knownCanonical_{tss_tes}_windows_201_100.bed",
        comb="results/XR/{samples}/{samples}_hg19_sorted_chr.bed",
    output:
        TS=temp("results/XR/{samples}/{samples}_hg19_sorted_TS_{tss_tes}_windows.bed"),
        NTS=temp("results/XR/{samples}/{samples}_hg19_sorted_NTS_{tss_tes}_windows.bed"),
        tss=temp("results/XR/{samples}/{samples}_hg19_sorted_{tss_tes}.bed"),
        tss_comb=temp("results/XR/{samples}/{samples}_hg19_sorted_{tss_tes}_combined.bed"),
        info=temp("results/XR/{samples}/{samples}_hg19_sorted_{tss_tes}_combined_full.bed"),
        rpkm="results/XR/{samples}/{samples}_hg19_sorted_{tss_tes}_combined_rpkm.bed",
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[1]),
    log:
        "logs/rule/analysis/{samples}/{samples}_hg19_{tss_tes}_xr.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_hg19_{tss_tes}_xr.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:
        """
        (echo "`date -R`: Counting reads on TS..." &&
        bedtools intersect \
        -sorted -a {input.genes} \
        -b {input.comb} \
        -wa -c -S -F 0.5 |
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t""TS""\\t"$7}}' \
        > {output.TS} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Counting reads on NTS..." &&
        bedtools intersect \
        -sorted -a {input.genes} \
        -b {input.comb} \
        -wa -c -s -F 0.5 |
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t""NTS""\\t"$7}}' \
        > {output.NTS} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Cat TS & NTS files..." &&
        cat {output.TS} {output.NTS} > {output.tss} &&
        echo "`date -R`: Success! TS & NTS files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        
        (echo "`date -R`: Combine windows..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {output.tss} \
        -strand T \
        -o {output.tss_comb} &&
        echo "`date -R`: Success! Windows are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        moreinfo="$(echo {params.info} | sed 's/,/\\t/g' )"

        (echo "`date -R`: Add info..." &&
        python3 workflow/scripts/addColumns.py \
        -i {output.tss_comb} \
        -o {output.info} \
        -c $moreinfo "." {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Calculating RPKM values..." &&
        python3 workflow/scripts/RPKM.py \
        -i {output.info} \
        -chse 2 3 \
        -c 7 \
        -mr 0 \
        -o {output.rpkm} &&
        echo "`date -R`: Success! RPKMs are calculated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule tss_ds:
    input:
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes_knownCanonical_{tss_tes}_windows_201_100.bed",
        comb="results/DS/{samples}/{samples}_hg19_sorted_dipyrimidines_chr.bed",
    output:
        TS=temp("results/DS/{samples}/{samples}_hg19_sorted_dipyrimidines_TS_{tss_tes}_windows.bed"),
        NTS=temp("results/DS/{samples}/{samples}_hg19_sorted_dipyrimidines_NTS_{tss_tes}_windows.bed"),
        tss="results/DS/{samples}/{samples}_hg19_sorted_dipyrimidines_{tss_tes}.bed",
        tss_comb=temp("results/DS/{samples}/{samples}_hg19_sorted_dipyrimidines_{tss_tes}_combined.bed"),
        info=temp("results/DS/{samples}/{samples}_hg19_sorted_dipyrimidines_{tss_tes}_combined_full.bed"),
        rpkm="results/DS/{samples}/{samples}_hg19_sorted_dipyrimidines_{tss_tes}_combined_rpkm.bed",
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[1]),
    log:
        "logs/rule/analysis/{samples}/{samples}_hg19_{tss_tes}_ds.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_hg19_{tss_tes}_ds.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:
        """
        (echo "`date -R`: Counting reads on TS..." &&
        bedtools intersect \
        -sorted -a {input.genes} \
        -b {input.comb} \
        -wa -c -S -F 0.5 |
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t""TS""\\t"$7}}' \
        > {output.TS} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Counting reads on NTS..." &&
        bedtools intersect \
        -sorted -a {input.genes} \
        -b {input.comb} \
        -wa -c -s -F 0.5 |
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t""NTS""\\t"$7}}' \
        > {output.NTS} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Cat TS & NTS files..." &&
        cat {output.TS} {output.NTS} > {output.tss} &&
        echo "`date -R`: Success! TS & NTS files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
                
        (echo "`date -R`: Combine windows..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {output.tss} \
        -strand T \
        -o {output.tss_comb} &&
        echo "`date -R`: Success! Windows are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        moreinfo="$(echo {params.info} | sed 's/,/\\t/g' )"

        (echo "`date -R`: Add info..." &&
        python3 workflow/scripts/addColumns.py \
        -i {output.tss_comb} \
        -o {output.info} \
        -c $moreinfo "." {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Calculating RPKM values..." &&
        python3 workflow/scripts/RPKM.py \
        -i {output.info} \
        -chse 2 3 \
        -c 7 \
        -mr 0 \
        -o {output.rpkm} &&
        echo "`date -R`: Success! RPKMs are calculated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
rule tss_xr_sim:
    input:
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes_knownCanonical_{tss_tes}_windows_201_100.bed",
        comb="results/sim/{samples}/{samples}_hg19_sorted_chr.bed",
    output:
        TS=temp("results/sim/{samples}/{samples}_hg19_sorted_TS_{tss_tes}_windows.bed"),
        NTS=temp("results/sim/{samples}/{samples}_hg19_sorted_NTS_{tss_tes}_windows.bed"),
        tss=temp("results/sim/{samples}/{samples}_hg19_sorted_{tss_tes}.bed"),
        tss_comb=temp("results/sim/{samples}/{samples}_hg19_sorted_{tss_tes}_combined.bed"),
        info=temp("results/sim/{samples}/{samples}_hg19_sorted_{tss_tes}_combined_full.bed"),
        rpkm="results/sim/{samples}/{samples}_hg19_sorted_{tss_tes}_combined_rpkm.bed",
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[1]),
    log:
        "logs/rule/analysis/{samples}/{samples}_sim_hg19_{tss_tes}_xr.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_sim_hg19_{tss_tes}_xr.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:
        """
        (echo "`date -R`: Counting reads on TS..." &&
        bedtools intersect \
        -sorted -a {input.genes} \
        -b {input.comb} \
        -wa -c -S -F 0.5 |
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t""TS""\\t"$7}}' \
        > {output.TS} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Counting reads on NTS..." &&
        bedtools intersect \
        -sorted -a {input.genes} \
        -b {input.comb} \
        -wa -c -s -F 0.5 |
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t""NTS""\\t"$7}}' \
        > {output.NTS} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Cat TS & NTS files..." &&
        cat {output.TS} {output.NTS} > {output.tss} &&
        echo "`date -R`: Success! TS & NTS files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        
        (echo "`date -R`: Combine windows..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {output.tss} \
        -strand T \
        -o {output.tss_comb} &&
        echo "`date -R`: Success! Windows are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        moreinfo="$(echo {params.info} | sed 's/,/\\t/g' )"

        (echo "`date -R`: Add info..." &&
        python3 workflow/scripts/addColumns.py \
        -i {output.tss_comb} \
        -o {output.info} \
        -c $moreinfo "." {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Calculating RPKM values..." &&
        python3 workflow/scripts/RPKM.py \
        -i {output.info} \
        -chse 2 3 \
        -c 7 \
        -mr 0 \
        -o {output.rpkm} &&
        echo "`date -R`: Success! RPKMs are calculated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule tss_ds_sim:
    input:
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes_knownCanonical_{tss_tes}_windows_201_100.bed",
        comb="results/sim/{samples}/{samples}_hg19_sorted_dipyrimidines_chr.bed",
    output:
        TS=temp("results/sim/{samples}/{samples}_hg19_sorted_dipyrimidines_TS_{tss_tes}_windows.bed"),
        NTS=temp("results/sim/{samples}/{samples}_hg19_sorted_dipyrimidines_NTS_{tss_tes}_windows.bed"),
        tss="results/sim/{samples}/{samples}_hg19_sorted_dipyrimidines_{tss_tes}.bed",
        tss_comb=temp("results/sim/{samples}/{samples}_hg19_sorted_dipyrimidines_{tss_tes}_combined.bed"),
        info=temp("results/sim/{samples}/{samples}_hg19_sorted_dipyrimidines_{tss_tes}_combined_full.bed"),
        rpkm="results/sim/{samples}/{samples}_hg19_sorted_dipyrimidines_{tss_tes}_combined_rpkm.bed",
    params:
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[1]),
    log:
        "logs/rule/analysis/{samples}/{samples}_sim_hg19_{tss_tes}_ds.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_sim_hg19_{tss_tes}_ds.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:
        """
        (echo "`date -R`: Counting reads on TS..." &&
        bedtools intersect \
        -sorted -a {input.genes} \
        -b {input.comb} \
        -wa -c -S -F 0.5 |
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t""TS""\\t"$7}}' \
        > {output.TS} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Counting reads on NTS..." &&
        bedtools intersect \
        -sorted -a {input.genes} \
        -b {input.comb} \
        -wa -c -s -F 0.5 |
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t""NTS""\\t"$7}}' \
        > {output.NTS} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Cat TS & NTS files..." &&
        cat {output.TS} {output.NTS} > {output.tss} &&
        echo "`date -R`: Success! TS & NTS files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
                
        (echo "`date -R`: Combine windows..." &&
        python3 workflow/scripts/combinewindows.py \
        -i {output.tss} \
        -strand T \
        -o {output.tss_comb} &&
        echo "`date -R`: Success! Windows are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        moreinfo="$(echo {params.info} | sed 's/,/\\t/g' )"

        (echo "`date -R`: Add info..." &&
        python3 workflow/scripts/addColumns.py \
        -i {output.tss_comb} \
        -o {output.info} \
        -c $moreinfo "." {params.mappedReads} &&
        echo "`date -R`: Success! Info is added." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Calculating RPKM values..." &&
        python3 workflow/scripts/RPKM.py \
        -i {output.info} \
        -chse 2 3 \
        -c 7 \
        -mr 0 \
        -o {output.rpkm} &&
        echo "`date -R`: Success! RPKMs are calculated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
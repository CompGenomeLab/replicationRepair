rule pre_mapping_xr:
    input:
        plus="resources/samples/XR/{samples}_{build}_sorted_plus.bed",
        minus="resources/samples/XR/{samples}_{build}_sorted_minus.bed",
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes.bed",
    output:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_damSite.bed",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_damSite.bed",
        plus_intergenic="results/XR/{samples}/{samples}_{build}_sorted_xr_intergenic_plus.bed",
        minus_intergenic="results/XR/{samples}/{samples}_{build}_sorted_xr_intergenic_minus.bed",
    params:
        filt="'^chr([1-9]|1[0-9]|2[0-2]|X)'",
    log:
        "logs/{samples}/{samples}_{build}_pre_mapping_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_pre_mapping_xr.benchmark.txt",   
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """   
        (echo "`date -R`: Extracting exact damage sites (plus strand)..." &&
        awk '{{print $1"\\t"$3-8"\\t"$3-7"\\t""reads""\\t"".""\\t"$6}}' {input.plus} |
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.plus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Extracting exact damage sites (minus strand)..." &&
        awk '{{print $1"\\t"$2+7"\t"$2+8"\\t""reads""\\t"".""\\t"$6}}' {input.minus} |
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.minus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Getting intergenic (plus strand)..." &&
        bedtools intersect \
        -a {output.plus} \
        -b {input.genes} \
        -v -f 0.5 |&
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.plus_intergenic} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Getting intergenic (minus strand)..." &&
        bedtools intersect \
        -a {output.minus} \
        -b {input.genes} \
        -v -f 0.5 |&
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.minus_intergenic} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule pre_mapping_ds:
    input:
        plus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed", 
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes.bed",
    output:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_damSite.bed",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_damSite.bed",
        plus_intergenic="results/DS/{samples}/{samples}_{build}_sorted_ds_intergenic_plus.bed",
        minus_intergenic="results/DS/{samples}/{samples}_{build}_sorted_ds_intergenic_minus.bed",
    params:
        filt="'^chr([1-9]|1[0-9]|2[0-2]|X)'",  
    log:
        "logs/{samples}/{samples}_{build}_pre_mapping_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_pre_mapping_ds.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """   
        (echo "`date -R`: Extracting exact damage sites (plus strand)..." &&
        awk '{{print $1"\\t"$3-5"\t"$3-4"\\t""reads""\\t"".""\\t"$6}}' {input.plus} |
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.plus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Extracting exact damage sites (minus strand)..." &&
        awk '{{print $1"\\t"$2+4"\t"$2+5"\\t""reads""\\t"".""\\t"$6}}' {input.minus} |
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.minus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Getting intergenic (plus strand)..." &&
        bedtools intersect \
        -a {output.plus} \
        -b {input.genes} \
        -v -f 0.5 |&
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.plus_intergenic} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Getting intergenic (minus strand)..." &&
        bedtools intersect \
        -a {output.minus} \
        -b {input.genes} \
        -v -f 0.5 |&
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.minus_intergenic} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule pre_mapping_xr_sim:
    input:
        bed="resources/samples/sim/{samples}_{build}_xr_sim.bed", 
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes.bed",
    output:
        plus=temp("results/sim/{samples}/{samples}_{build}_xr_sim_plus_sorted.bed"),
        minus=temp("results/sim/{samples}/{samples}_{build}_xr_sim_minus_sorted.bed"),
        plus_dam="results/sim/{samples}/{samples}_{build}_xr_sim_plus_damSite.bed",
        minus_dam="results/sim/{samples}/{samples}_{build}_xr_sim_minus_damSite.bed",
        plus_intergenic="results/sim/{samples}/{samples}_{build}_sorted_xr_sim_intergenic_plus.bed",
        minus_intergenic="results/sim/{samples}/{samples}_{build}_sorted_xr_sim_intergenic_minus.bed",
    params:
        filt="'^chr([1-9]|1[0-9]|2[0-2]|X)'",    
    log:
        "logs/{samples}/{samples}_{build}_pre_mapping_sim_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_pre_mapping_sim_xr.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """   
        (echo "`date -R`: Separating plus strands..." &&
        awk '{{if($6=="+"){{print}}}}' {input.bed} > {output.plus} &&
        echo "`date -R`: Success! Strands are separated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Separating minus strands..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Strands are separated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Extracting exact damage sites (plus strand)..." &&
        awk '{{print $1"\\t"$3-8"\\t"$3-7"\\t""reads""\\t"".""\\t"$6}}' {output.plus} |
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.plus_dam} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus_dam}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Extracting exact damage sites (minus strand)..." &&
        awk '{{print $1"\\t"$2+7"\t"$2+8"\\t""reads""\\t"".""\\t"$6}}' {output.minus} |
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.minus_dam} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus_dam}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Getting intergenic (plus strand)..." &&
        bedtools intersect \
        -a {output.plus_dam} \
        -b {input.genes} \
        -v -f 0.5 |&
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.plus_intergenic} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Getting intergenic (minus strand)..." &&
        bedtools intersect \
        -a {output.minus_dam} \
        -b {input.genes} \
        -v -f 0.5 |&
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.minus_intergenic} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule pre_mapping_ds_sim:
    input:
        bed="resources/samples/sim/{samples}_{build}_ds_sim.bed",
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes.bed",
    output:
        plus=temp("results/sim/{samples}/{samples}_{build}_ds_sim_plus_sorted.bed"),
        minus=temp("results/sim/{samples}/{samples}_{build}_ds_sim_minus_sorted.bed"),
        plus_dam="results/sim/{samples}/{samples}_{build}_ds_sim_plus_damSite.bed",
        minus_dam="results/sim/{samples}/{samples}_{build}_ds_sim_minus_damSite.bed",
        plus_intergenic="results/sim/{samples}/{samples}_{build}_sorted_ds_sim_intergenic_plus.bed",
        minus_intergenic="results/sim/{samples}/{samples}_{build}_sorted_ds_sim_intergenic_minus.bed",
    params:
        filt="'^chr([1-9]|1[0-9]|2[0-2]|X)'",    
    log:
        "logs/{samples}/{samples}_{build}_pre_mapping_sim_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_pre_mapping_sim_ds.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """   
        (echo "`date -R`: Separating plus strands..." &&
        awk '{{if($6=="+"){{print}}}}' {input.bed} > {output.plus} &&
        echo "`date -R`: Success! Strands are separated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Separating minus strands..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Strands are separated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Extracting exact damage sites (plus strand)..." &&
        awk '{{print $1"\\t"$3-5"\t"$3-4"\\t""reads""\\t"".""\\t"$6}}' {output.plus} |
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.plus_dam} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus_dam}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Extracting exact damage sites (minus strand)..." &&
        awk '{{print $1"\\t"$2+4"\t"$2+5"\\t""reads""\\t"".""\\t"$6}}' {output.minus} |
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.minus_dam} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus_dam}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Getting intergenic (plus strand)..." &&
        bedtools intersect \
        -a {output.plus_dam} \
        -b {input.genes} \
        -v -f 0.5 |&
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.plus_intergenic} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Getting intergenic (minus strand)..." &&
        bedtools intersect \
        -a {output.minus_dam} \
        -b {input.genes} \
        -v -f 0.5 |&
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.minus_intergenic} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
       """

rule pre_mapping_chipseq:
    input:
        bed=lambda w: input4chip(w, config["chipseq"]["samples"], 
            config["chipseq"]["srr"]["enabled"], config["chipseq"]["srr"]["codes"]),
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes.bed",
    output:
        org="results/chipseq/{samples}/{samples}_{build}_org.bed",
        intergenic="results/chipseq/{samples}/{samples}_{build}_org_intergenic.bed",
    params:
        filt="'^chr([1-9]|1[0-9]|2[0-2]|X)'",  
        marker_name=lambda w: getMarkerName(w), 
    log:
        "logs/{samples}/{samples}_{build}_pre_mapping_chipseq.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_pre_mapping_chipseq.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Reorganizing..." &&
        awk '{{print $1"\\t"$2"\\t"$3"\\t""{params.marker_name}""\\t"".""\\t"$6}}' {input.bed} |
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.org} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.org}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Getting intergenic..." &&
        bedtools intersect \
        -a {output.org} \
        -b {input.genes} \
        -v -f 0.5 |&
        sort -u -k1,1 -k2,2n -k3,3n |&
        egrep {params.filt} > {output.intergenic} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule pre_mapping_methyl:
    input:
        bed="resources/samples/methyl/{samples}_{build}_chr.bed",
        genes="resources/ref_genomes/hg19/hg19_ucsc_genes.bed",
    output:
        intergenic="results/methyl/{samples}/{samples}_{build}_org_intergenic.bed",
    params:
        filt="'^chr([1-9]|1[0-9]|2[0-2]|X)'",  
        marker_name=lambda w: getMarkerName(w, idx=5), 
    log:
        "logs/{samples}/{samples}_{build}_pre_mapping_methyl.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_pre_mapping_methyl.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Getting intergenic..." &&
        bedtools intersect \
        -a {input.bed} \
        -b {input.genes} \
        -v -f 0.5 \
        > {output.intergenic} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """
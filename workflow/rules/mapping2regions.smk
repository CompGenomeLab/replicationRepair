#!/bin/env python

rule mapping2regions_xr:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_damSite.bed",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_damSite.bed",
    output:
        sorted_region=temp("results/XR/{samples}/{samples}_{build}_{regions}_sorted.bed"),
        int_plus=temp("results/XR/{samples}/{samples}_{build}_plus_{regions}.txt"),
        int_minus=temp("results/XR/{samples}/{samples}_{build}_minus_{regions}.txt"),
        comb_plus=temp("results/XR/{samples}/{samples}_{build}_plus_{regions}_combined.txt"),
        comb_minus=temp("results/XR/{samples}/{samples}_{build}_minus_{regions}_combined.txt"),
        info_plus=temp("results/XR/{samples}/{samples}_{build}_plus_{regions}_combined_info.txt"),
        info_minus=temp("results/XR/{samples}/{samples}_{build}_minus_{regions}_combined_info.txt"),
        rpkm_plus="results/XR/{samples}/{samples}_{build}_plus_{regions}_combined_rpkm.txt",
        rpkm_minus="results/XR/{samples}/{samples}_{build}_minus_{regions}_combined_rpkm.txt",
    params:
        sdir="results/XR/{samples}",
        sname="{samples}_{build}",
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
        region_name="{regions}",
        getComb=lambda w: getCombine(w.regions, config["region_comb_opt"], config["regions"]), 
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[0], input[1]),
    log:
        "logs/XR/{samples}/{samples}_{build}_mapping2regions_xr_{regions}.log",
    benchmark:
        "logs/XR/{samples}/{samples}_{build}_mapping2regions_xr_{regions}.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        workflow/scripts/mapping2regions.sh \
        {input.plus} \
        {input.minus} \
        {params.sdir} \
        {params.sname} \
        {params.region} \
        {params.region_name} \
        {params.getComb} \
        {params.info} \
        {params.mappedReads} \
        {log}
        """

rule mapping2regions_ds:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_damSite.bed",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_damSite.bed",
    output:
        sorted_region=temp("results/DS/{samples}/{samples}_{build}_{regions}_sorted.bed"),
        int_plus=temp("results/DS/{samples}/{samples}_{build}_plus_{regions}.txt"),
        int_minus=temp("results/DS/{samples}/{samples}_{build}_minus_{regions}.txt"),
        comb_plus=temp("results/DS/{samples}/{samples}_{build}_plus_{regions}_combined.txt"),
        comb_minus=temp("results/DS/{samples}/{samples}_{build}_minus_{regions}_combined.txt"),
        info_plus=temp("results/DS/{samples}/{samples}_{build}_plus_{regions}_combined_info.txt"),
        info_minus=temp("results/DS/{samples}/{samples}_{build}_minus_{regions}_combined_info.txt"),
        rpkm_plus="results/DS/{samples}/{samples}_{build}_plus_{regions}_combined_rpkm.txt",
        rpkm_minus="results/DS/{samples}/{samples}_{build}_minus_{regions}_combined_rpkm.txt",
    params:
        sdir="results/DS/{samples}",
        sname="{samples}_{build}",
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
        region_name="{regions}",
        getComb=lambda w: getCombine(w.regions, config["region_comb_opt"], config["regions"]), 
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[0], input[1]),
    log:
        "logs/DS/{samples}/{samples}_{build}_mapping2regions_ds_{regions}.log",
    benchmark:
        "logs/DS/{samples}/{samples}_{build}_mapping2regions_ds_{regions}.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        workflow/scripts/mapping2regions.sh \
        {input.plus} \
        {input.minus} \
        {params.sdir} \
        {params.sname} \
        {params.region} \
        {params.region_name} \
        {params.getComb} \
        {params.info} \
        {params.mappedReads} \
        {log}
        """

rule mapping2regions_xr_intergenic:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_intergenic_plus.bed",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_intergenic_minus.bed",
    output:
        sorted_region=temp("results/intergenic/{samples}/{samples}_{build}_{regions}_sorted.bed"),
        int_plus=temp("results/intergenic/{samples}/{samples}_{build}_plus_{regions}.txt"),
        int_minus=temp("results/intergenic/{samples}/{samples}_{build}_minus_{regions}.txt"),
        comb_plus=temp("results/intergenic/{samples}/{samples}_{build}_plus_{regions}_combined.txt"),
        comb_minus=temp("results/intergenic/{samples}/{samples}_{build}_minus_{regions}_combined.txt"),
        info_plus=temp("results/intergenic/{samples}/{samples}_{build}_plus_{regions}_combined_info.txt"),
        info_minus=temp("results/intergenic/{samples}/{samples}_{build}_minus_{regions}_combined_info.txt"),
        rpkm_plus="results/intergenic/{samples}/{samples}_{build}_plus_{regions}_combined_rpkm.txt",
        rpkm_minus="results/intergenic/{samples}/{samples}_{build}_minus_{regions}_combined_rpkm.txt",
    params:
        sdir="results/intergenic/{samples}",
        sname="{samples}_{build}",
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
        region_name="{regions}",
        getComb=lambda w: getCombine(w.regions, config["region_comb_opt"], config["regions"]), 
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[0], input[1]),
    log:
        "logs/XR/{samples}/{samples}_{build}_mapping2regions_xr_intergenic_{regions}.log",
    benchmark:
        "logs/XR/{samples}/{samples}_{build}_mapping2regions_xr_intergenic_{regions}.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        workflow/scripts/mapping2regions.sh \
        {input.plus} \
        {input.minus} \
        {params.sdir} \
        {params.sname} \
        {params.region} \
        {params.region_name} \
        {params.getComb} \
        {params.info} \
        {params.mappedReads} \
        {log}
        """

rule mapping2regions_ds_intergenic:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_intergenic_plus.bed",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_intergenic_minus.bed",
    output:
        sorted_region=temp("results/intergenic/{samples}/{samples}_{build}_{regions}_sorted.bed"),
        int_plus=temp("results/intergenic/{samples}/{samples}_{build}_plus_{regions}.txt"),
        int_minus=temp("results/intergenic/{samples}/{samples}_{build}_minus_{regions}.txt"),
        comb_plus=temp("results/intergenic/{samples}/{samples}_{build}_plus_{regions}_combined.txt"),
        comb_minus=temp("results/intergenic/{samples}/{samples}_{build}_minus_{regions}_combined.txt"),
        info_plus=temp("results/intergenic/{samples}/{samples}_{build}_plus_{regions}_combined_info.txt"),
        info_minus=temp("results/intergenic/{samples}/{samples}_{build}_minus_{regions}_combined_info.txt"),
        rpkm_plus="results/intergenic/{samples}/{samples}_{build}_plus_{regions}_combined_rpkm.txt",
        rpkm_minus="results/intergenic/{samples}/{samples}_{build}_minus_{regions}_combined_rpkm.txt",
    params:
        sdir="results/intergenic/{samples}",
        sname="{samples}_{build}",
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
        region_name="{regions}",
        getComb=lambda w: getCombine(w.regions, config["region_comb_opt"], config["regions"]), 
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[0], input[1]),
    log:
        "logs/DS/{samples}/{samples}_{build}_mapping2regions_ds_intergenic_{regions}.log",
    benchmark:
        "logs/DS/{samples}/{samples}_{build}_mapping2regions_ds_intergenic_{regions}.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        workflow/scripts/mapping2regions.sh \
        {input.plus} \
        {input.minus} \
        {params.sdir} \
        {params.sname} \
        {params.region} \
        {params.region_name} \
        {params.getComb} \
        {params.info} \
        {params.mappedReads} \
        {log}
        """

rule mapping2regions_xr_sim:
    input:
        plus="results/sim/{samples}/{samples}_{build}_xr_sim_plus_damSite.bed",
        minus="results/sim/{samples}/{samples}_{build}_xr_sim_minus_damSite.bed",
    output:
        sorted_region=temp("results/sim/{samples}/{samples}_{build}_{regions}_sorted.bed"),
        int_plus=temp("results/sim/{samples}/{samples}_{build}_plus_{regions}.txt"),
        int_minus=temp("results/sim/{samples}/{samples}_{build}_minus_{regions}.txt"),
        comb_plus=temp("results/sim/{samples}/{samples}_{build}_plus_{regions}_combined.txt"),
        comb_minus=temp("results/sim/{samples}/{samples}_{build}_minus_{regions}_combined.txt"),
        info_plus=temp("results/sim/{samples}/{samples}_{build}_plus_{regions}_combined_info.txt"),
        info_minus=temp("results/sim/{samples}/{samples}_{build}_minus_{regions}_combined_info.txt"),
        rpkm_plus="results/sim/{samples}/{samples}_{build}_plus_{regions}_combined_rpkm.txt",
        rpkm_minus="results/sim/{samples}/{samples}_{build}_minus_{regions}_combined_rpkm.txt",
    params:
        sdir="results/sim/{samples}",
        sname="{samples}_{build}",
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
        region_name="{regions}",
        getComb=lambda w: getCombine(w.regions, config["region_comb_opt"], config["regions"]), 
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[0], input[1]),
    log:
        "logs/XR/{samples}/{samples}_{build}_sim_mapping2regions_xr_{regions}.log",
    benchmark:
        "logs/XR/{samples}/{samples}_{build}_sim_mapping2regions_xr_{regions}.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        workflow/scripts/mapping2regions.sh \
        {input.plus} \
        {input.minus} \
        {params.sdir} \
        {params.sname} \
        {params.region} \
        {params.region_name} \
        {params.getComb} \
        {params.info} \
        {params.mappedReads} \
        {log}
        """

rule mapping2regions_ds_sim:
    input:
        plus="results/sim/{samples}/{samples}_{build}_ds_sim_plus_damSite.bed",
        minus="results/sim/{samples}/{samples}_{build}_ds_sim_minus_damSite.bed",
    output:
        sorted_region=temp("results/sim/{samples}/{samples}_{build}_{regions}_sorted.bed"),
        int_plus=temp("results/sim/{samples}/{samples}_{build}_plus_{regions}.txt"),
        int_minus=temp("results/sim/{samples}/{samples}_{build}_minus_{regions}.txt"),
        comb_plus=temp("results/sim/{samples}/{samples}_{build}_plus_{regions}_combined.txt"),
        comb_minus=temp("results/sim/{samples}/{samples}_{build}_minus_{regions}_combined.txt"),
        info_plus=temp("results/sim/{samples}/{samples}_{build}_plus_{regions}_combined_info.txt"),
        info_minus=temp("results/sim/{samples}/{samples}_{build}_minus_{regions}_combined_info.txt"),
        rpkm_plus="results/sim/{samples}/{samples}_{build}_plus_{regions}_combined_rpkm.txt",
        rpkm_minus="results/sim/{samples}/{samples}_{build}_minus_{regions}_combined_rpkm.txt",
    params:
        sdir="results/sim/{samples}",
        sname="{samples}_{build}",
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
        region_name="{regions}",
        getComb=lambda w: getCombine(w.regions, config["region_comb_opt"], config["regions"]),
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[0], input[1]), 
    log:
        "logs/DS/{samples}/{samples}_{build}_sim_mapping2regions_ds_{regions}.log",
    benchmark:
        "logs/DS/{samples}/{samples}_{build}_sim_mapping2regions_ds_{regions}.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        workflow/scripts/mapping2regions.sh \
        {input.plus} \
        {input.minus} \
        {params.sdir} \
        {params.sname} \
        {params.region} \
        {params.region_name} \
        {params.getComb} \
        {params.info} \
        {params.mappedReads} \
        {log}
        """

rule mapping2regions_xr_intergenic_sim:
    input:
        plus="results/sim/{samples}/{samples}_{build}_sorted_xr_sim_intergenic_plus.bed",
        minus="results/sim/{samples}/{samples}_{build}_sorted_xr_sim_intergenic_minus.bed",
    output:
        sorted_region=temp("results/intergenic_sim/{samples}/{samples}_{build}_{regions}_sorted.bed"),
        int_plus=temp("results/intergenic_sim/{samples}/{samples}_{build}_plus_{regions}.txt"),
        int_minus=temp("results/intergenic_sim/{samples}/{samples}_{build}_minus_{regions}.txt"),
        comb_plus=temp("results/intergenic_sim/{samples}/{samples}_{build}_plus_{regions}_combined.txt"),
        comb_minus=temp("results/intergenic_sim/{samples}/{samples}_{build}_minus_{regions}_combined.txt"),
        info_plus=temp("results/intergenic_sim/{samples}/{samples}_{build}_plus_{regions}_combined_info.txt"),
        info_minus=temp("results/intergenic_sim/{samples}/{samples}_{build}_minus_{regions}_combined_info.txt"),
        rpkm_plus="results/intergenic_sim/{samples}/{samples}_{build}_plus_{regions}_combined_rpkm.txt",
        rpkm_minus="results/intergenic_sim/{samples}/{samples}_{build}_minus_{regions}_combined_rpkm.txt",
    params:
        sdir="results/intergenic_sim/{samples}",
        sname="{samples}_{build}",
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
        region_name="{regions}",
        getComb=lambda w: getCombine(w.regions, config["region_comb_opt"], config["regions"]), 
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[0], input[1]),
    log:
        "logs/XR/{samples}/{samples}_{build}_sim_mapping2regions_xr_intergenic_{regions}.log",
    benchmark:
        "logs/XR/{samples}/{samples}_{build}_sim_mapping2regions_xr_intergenic_{regions}.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        workflow/scripts/mapping2regions.sh \
        {input.plus} \
        {input.minus} \
        {params.sdir} \
        {params.sname} \
        {params.region} \
        {params.region_name} \
        {params.getComb} \
        {params.info} \
        {params.mappedReads} \
        {log}
        """

rule mapping2regions_ds_intergenic_sim:
    input:
        plus="results/sim/{samples}/{samples}_{build}_sorted_ds_sim_intergenic_plus.bed",
        minus="results/sim/{samples}/{samples}_{build}_sorted_ds_sim_intergenic_minus.bed",
    output:
        sorted_region=temp("results/intergenic_sim/{samples}/{samples}_{build}_{regions}_sorted.bed"),
        int_plus=temp("results/intergenic_sim/{samples}/{samples}_{build}_plus_{regions}.txt"),
        int_minus=temp("results/intergenic_sim/{samples}/{samples}_{build}_minus_{regions}.txt"),
        comb_plus=temp("results/intergenic_sim/{samples}/{samples}_{build}_plus_{regions}_combined.txt"),
        comb_minus=temp("results/intergenic_sim/{samples}/{samples}_{build}_minus_{regions}_combined.txt"),
        info_plus=temp("results/intergenic_sim/{samples}/{samples}_{build}_plus_{regions}_combined_info.txt"),
        info_minus=temp("results/intergenic_sim/{samples}/{samples}_{build}_minus_{regions}_combined_info.txt"),
        rpkm_plus="results/intergenic_sim/{samples}/{samples}_{build}_plus_{regions}_combined_rpkm.txt",
        rpkm_minus="results/intergenic_sim/{samples}/{samples}_{build}_minus_{regions}_combined_rpkm.txt",
    params:
        sdir="results/intergenic_sim/{samples}",
        sname="{samples}_{build}",
        region=lambda w: getRegion(w.regions, config["region_file"], config["regions"]),
        region_name="{regions}",
        getComb=lambda w: getCombine(w.regions, config["region_comb_opt"], config["regions"]),
        info=lambda w: info(w),
        mappedReads=lambda w, input: mappedReads(input[0], input[1]), 
    log:
        "logs/DS/{samples}/{samples}_{build}_sim_mapping2regions_ds_intergenic_{regions}.log",
    benchmark:
        "logs/DS/{samples}/{samples}_{build}_sim_mapping2regions_ds_intergenic_{regions}.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        workflow/scripts/mapping2regions.sh \
        {input.plus} \
        {input.minus} \
        {params.sdir} \
        {params.sname} \
        {params.region} \
        {params.region_name} \
        {params.getComb} \
        {params.info} \
        {params.mappedReads} \
        {log}
        """
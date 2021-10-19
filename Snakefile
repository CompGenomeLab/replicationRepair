#!/bin/env python

#### configuration file #### 
configfile: "config/config.yaml"

include: "workflow/rules/common.smk"

wildcard_constraints:
    regions='|'.join([r for r in (config["regions"] + config["regions_mut"])]),
    build=config["build"],
    samples='|'.join([s for s in (config["sample_ds"] + config["sample_xr"] + config["sample_edu"] + config["sample_okseq"] + config["sample_mutation"] + config["chipseq"]["samples"])]),
    tss_tes='tss|tes'

rule all:
    input:
        #lambda w: allInput(config["build"], config["sample_edu"], "edu"),
        #lambda w: allInput(config["build"], config["sample_okseq"], "okseq"),
        #lambda w: allInput(config["build"], config["sample_mutation"], "mutation", config["regions_mut"]),
        lambda w: allInput(config["build"], config["chipseq"]["samples"], "markers_intergenic", config["regions"]),
        lambda w: allInput(config["build"], config["sample_ds"], "ds", config["regions"]),
        lambda w: allInput(config["build"], config["sample_xr"], "xr", config["regions"]),
        lambda w: allInput(build=config["build"], method="report", 
            regions=config["regions"]),

# prepare genome
include: "workflow/rules/genome_build.smk"
include: "workflow/rules/genome_indexing.smk"

# EdU
include: "workflow/snakefiles/EdU"

# OK-seq
#include: "workflow/snakefiles/OK-seq"

# chip-seq
if config["chipseq"]["srr"]["enabled"]:
    include: "workflow/rules/sra_chipseq.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/bam2bed.smk" 

# Mutation Analyses
include: "workflow/rules/get_sbs_muts.smk"
include: "workflow/rules/organize.smk"
include: "workflow/rules/filter_target_muts.smk"
include: "workflow/rules/sep_strands.smk"
include: "workflow/rules/intergenic.smk"

# Process Retrieved Data
#include: "workflow/rules/make_windows.smk"
#include: "workflow/rules/????.smk"

# Further Analyses
#include: "workflow/rules/sep_strands.smk"
include: "workflow/rules/ts_nts.smk"
include: "workflow/rules/tss.smk"
include: "workflow/rules/pre_mapping.smk"
include: "workflow/rules/mapping2regions.smk"
include: "workflow/rules/combine_files.smk"

# Plots
include: "workflow/rules/figure1.smk"
include: "workflow/rules/figure2.smk"
include: "workflow/rules/figure3.smk"
#!/bin/env python

#### configuration file #### 
configfile: "config/config.yaml"

include: "workflow/rules/common.smk"

wildcard_constraints:
    regions='|'.join([r for r in (config["regions"] + config["regions_mut"])]),
    build=config["build"],
    samples='|'.join([s for s in (config["sample_ds"] + config["sample_xr"] + config["sample_edu"] + config["sample_input"] + config["sample_okseq"] + config["sample_mutation"])]),
    tss_tes='tss|tes'

rule all:
    input:
        #lambda w: allInput(config["build"], config["sample_input"], "input"),
        lambda w: allInput(config["build"], config["sample_edu"], "edu"),
        #lambda w: allInput(config["build"], config["sample_okseq"], "okseq"),
        #lambda w: allInput(config["build"], config["sample_mutation"], "mutation", config["regions_mut"]),
        #lambda w: allInput(config["build"], config["sample_ds"], "ds", config["regions"]),
        #lambda w: allInput(config["build"], config["sample_xr"], "xr", config["regions"]),
        #lambda w: allInput(build=config["build"], method="report", 
        #    regions=config["regions"]),

# prepare genome
include: "workflow/rules/genome_build.smk"
include: "workflow/rules/genome_indexing.smk"
include: "workflow/rules/genome_idx2ron.smk"

# DNA-seq
#include: "workflow/snakefiles/DNA-seq"

# EdU
include: "workflow/snakefiles/EdU"

# OK-seq
#include: "workflow/snakefiles/OK-seq"

# Mutation Analyses
#include: "workflow/snakefiles/mutation"

# Process Retrieved Data
#include: "workflow/rules/make_windows.smk"
#include: "workflow/rules/????.smk"

# Further Analyses
#include: "workflow/rules/sep_strands.smk"
include: "workflow/rules/ts_nts.smk"
include: "workflow/rules/tss.smk"
include: "workflow/rules/pre_mapping.smk"
include: "workflow/rules/mapping2regions.smk"
#include: "workflow/rules/intergenic.smk"
#include: "workflow/rules/intersect.smk"
#include: "workflow/rules/combine_windows.smk"
#include: "workflow/rules/more_info.smk"
#include: "workflow/rules/rpkm.smk"
include: "workflow/rules/combine_files.smk"


# Plots
#include: "workflow/rules/????.smk"
#include: "workflow/rules/????.smk"
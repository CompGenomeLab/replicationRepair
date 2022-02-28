#!/bin/env python

#### configuration file #### 
configfile: "config/config.yaml"

include: "workflow/rules/common.smk"

wildcard_constraints:
    regions='|'.join([r for r in (config["regions"] + config["regions_mut"])]),
    build=config["build"],
    samples='|'.join([s for s in (
                                config["xr"]["samples"] + 
                                config["ds"]["samples"] + 
                                config["edu"]["samples"] + 
                                config["okseq"]["samples"] + 
                                config["mutation"]["samples"] 
                                )]),
    tss_tes='tss|tes'

rule all:
    input:
        lambda w: allInput(
            build=config["build"], 
            sampleList=config["edu"]["samples"], 
            srrEnabled=config["edu"]["srr"]["enabled"], 
            srrList=config["edu"]["srr"]["codes"],
            method="edu"
            ),
        lambda w: allInput(
            build=config["build"], 
            sampleList=config["okseq"]["samples"],
            srrEnabled=config["okseq"]["srr"]["enabled"], 
            srrList=config["okseq"]["srr"]["codes"], 
            method="okseq"
            ),
        lambda w: allInput(
            build=config["build"], 
            sampleList=config["mutation"]["samples"], 
            method="mutation", 
            regions=config["regions_mut"]
            ),
        lambda w: allInput(
            build=config["build"], 
            sampleList=config["ds"]["samples"],  
            method="ds", 
            regions=config["regions"]
            ),
        lambda w: allInput(
            build=config["build"], 
            sampleList=config["xr"]["samples"],  
            method="xr", 
            regions=config["regions"]
            ),
        lambda w: allInput(
            build=config["build"], 
            method="report", 
            regions=config["regions"]
            ),

# Prepare genome
#include: "workflow/rules/genome_build.smk"
#include: "workflow/rules/genome_indexing.smk"

# EdU Analysis
if config["edu"]["srr"]["enabled"]:
    include: "workflow/rules/sra_edu.smk"
include: "workflow/rules/fastqc.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/bam2bed.smk"
include: "workflow/rules/sort_filter.smk"
include: "workflow/rules/sep_strands.smk"
include: "workflow/rules/genomecov.smk"
include: "workflow/rules/bedGraphToBigWig.smk"
include: "workflow/rules/bam_correlation.smk"
include: "workflow/rules/replication_timing.smk"
include: "workflow/rules/produceReplicationDomains.smk"
include: "workflow/rules/getPublicRepli.smk"

# OK-seq Analysis
if config["okseq"]["srr"]["enabled"]:
    include: "workflow/rules/sra_okseq.smk"
#include: "workflow/rules/fastqc.smk" # The rule implemented in EdU Analysis
include: "workflow/rules/adaptor_handling.smk"
#include: "workflow/rules/align.smk" # The rule implemented in EdU Analysis
#include: "workflow/rules/bam_correlation.smk" # The rule implemented in EdU Analysis
include: "workflow/rules/produceInitiationZones.smk"
include: "workflow/rules/okseqIntersect.smk"
include: "workflow/rules/getHeLaIZ.smk"

# Mutation Analysis
include: "workflow/rules/get_sbs_muts.smk"
include: "workflow/rules/organize.smk"
include: "workflow/rules/filter_target_muts.smk"
#include: "workflow/rules/sep_strands.smk" # The rule implemented in EdU Analysis
include: "workflow/rules/intergenic.smk"

# Chromatin States
include: "workflow/rules/getChromHMM.smk"

# Process Retrieved Data
include: "workflow/rules/intersect2repDomains.smk"
include: "workflow/rules/make_windows.smk"
include: "workflow/rules/countMotifs.smk"

# Further Analyses
include: "workflow/rules/combine_replicates.smk"
include: "workflow/rules/ts_nts.smk"
include: "workflow/rules/tss.smk"
include: "workflow/rules/pre_mapping.smk"
include: "workflow/rules/mapping2regions.smk"
include: "workflow/rules/combine_files.smk"

# Figures
include: "workflow/rules/figure1.smk"
include: "workflow/rules/figure2.smk"
include: "workflow/rules/figure3_S5.smk"
include: "workflow/rules/figure4B.smk"
include: "workflow/rules/figure4C_4D_S10_S11_S12B_S12C.smk"
include: "workflow/rules/figure4E_S13_S14_S15_S16_S17.smk"
include: "workflow/rules/figure5A.smk"
include: "workflow/rules/figure5B_5C_5D.smk"
include: "workflow/rules/bam_corr_graphs.smk" # contains figure S2A
include: "workflow/rules/figureS3B_S3C_S3D.smk"
include: "workflow/rules/figureS4.smk"
include: "workflow/rules/figureS6_S7.smk"
include: "workflow/rules/figureS8_S9.smk"
include: "workflow/rules/figureS18.smk"
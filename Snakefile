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
                                config["mutation"]["samples"] + 
                                config["chipseq"]["samples"] +
                                config["methyl"]["samples"]
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
        #lambda w: allInput(
        #    build=config["build"], 
        #    sampleList=config["chipseq"]["samples"], 
        #    srrEnabled=config["chipseq"]["srr"]["enabled"], 
        #    srrList=config["chipseq"]["srr"]["codes"], 
        #    method="markers_intergenic", 
        #    regions=config["regions"]
        #    ),
        #lambda w: allInput(
        #    build=config["build"], 
        #    sampleList=config["methyl"]["samples"],  
        #    method="methyl", 
        #    regions=config["regions"]
        #    ),
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
include: "workflow/rules/genome_build.smk"
include: "workflow/rules/genome_indexing.smk"

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

# OK-seq Analysis
if config["okseq"]["srr"]["enabled"]:
    include: "workflow/rules/sra_okseq.smk"
include: "workflow/rules/fastqc.smk"
include: "workflow/rules/adaptor_handling.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/bam_correlation.smk"
include: "workflow/rules/produceInitiationZones.smk"
include: "workflow/rules/okseqIntersect.smk"
include: "workflow/rules/getHeLaIZ.smk"

# Chip-seq Analysis
if config["chipseq"]["srr"]["enabled"]:
    include: "workflow/rules/sra_chipseq.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/bam2bed.smk" 
include: "workflow/rules/bam_correlation.smk"

# Methylation
#include: "workflow/rules/liftMethyl.smk"

# Mutation Analysis
include: "workflow/rules/get_sbs_muts.smk"
include: "workflow/rules/organize.smk"
include: "workflow/rules/filter_target_muts.smk"
include: "workflow/rules/sep_strands.smk"
include: "workflow/rules/intergenic.smk"

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

# Plots
include: "workflow/rules/bam_corr_graphs.smk"
include: "workflow/rules/figure_markers.smk"
include: "workflow/rules/figure1.smk"
include: "workflow/rules/figure2.smk"
include: "workflow/rules/figure3.smk"
include: "workflow/rules/figure4_5.smk"
include: "workflow/rules/figure6A.smk"
include: "workflow/rules/figure6.smk"
include: "workflow/rules/figureS2.smk"
include: "workflow/rules/figureS3.smk"
include: "workflow/rules/figureS5.smk"
include: "workflow/rules/figureS5_repdomains.smk"
include: "workflow/rules/figureS8_S9.smk"
include: "workflow/rules/figure_seq_context_asymmetry.smk"
include: "workflow/rules/figure_normDSXR.smk"
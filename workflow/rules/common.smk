#!/bin/env python 

import warnings
import os
import subprocess

################### Helper Functions ###########################################

def isSingle(sample, sampleList, method):

    sample_dir = "resources/samples/" + method + "/"
    single = sample_dir + sample + ".fastq.gz"
    pairedR1 = sample_dir + sample + "_R1.fastq.gz"
    
    if os.path.isfile(single) and "_R1" not in single:
        return True
    elif os.path.isfile(pairedR1):
        return False
    else:
        raise(ValueError(sample + ": Sample not found..."))

def input4filter(wildcards, sampleList, method):

    if isSingle(wildcards.samples, sampleList, method):
        return ("results/" + method + "/{samples}/{samples}_{build}_se.bed")
    else:    
        return ("results/" + method + "/{samples}/{samples}_{build}_pe.bed")

def input4peakCalling(wildcards, build, method):

    samp = "results/" + method + "/" + wildcards.samples + "/" + wildcards.samples + "_" + build + "_sorted.bam"

    if "rep2" in wildcards.samples:
        inp = samp.replace("rep2-", "rep2-inp")  
    else:
        inp = samp.replace("rep-", "rep-inp")  

    return [samp, inp]


def lineNum(file):
    
    linenum = 0
    if os.path.exists(file):
        with open(file) as f:
            for line in f:
                linenum += 1

    warnMessage = ("\n" + file + " file is either empty or does not exists!\n" + 
        "It is expected if this is a dry-run. The file will be produced " + 
        "after the execution.")

    if linenum == 0:
        warnings.warn(warnMessage)

    return linenum

def mappedReads(*files):

    lineNumber = 0
    for file in files:
        lineNumber += lineNum(str(file))

    return lineNumber

def info(wildcards, method="sample"):

    sample_name = wildcards.samples.replace("-","_")

    if method == "sample":

        with open("resources/samples/samples.csv") as f:

            targetLine = ""
            for line in f:

                if line.strip().split(",")[1] == sample_name:

                    #targetLine = line.strip().replace(",","\t")
                    targetLine = line.strip()
                    return targetLine
                
            if targetLine == "":
                raise(ValueError(sample_name + " not found in samples.csv file..."))
    
    elif method == "marker":

        with open("resources/samples/markers.csv") as f:

            targetLine = ""
            for line in f:

                if line.strip().split(",")[2] == sample_name:

                    #targetLine = line.strip().replace(",","\t")
                    targetLine = line.strip()
                    return targetLine
                
            if targetLine == "":
                raise(ValueError(sample_name + " not found in markers.csv file..."))

def getMarkerName(wildcards):
    return wildcards.samples.split("_")[1]


def getCombine(region, combineList, regionList):

    try:
        idx = regionList.index(region)
    except:
       raise(ValueError("Designated wildcard cannot be found in region list."))
        
    return combineList[idx]

def getRegion(region, rawRegionList, regionList):

    try:
        idx = regionList.index(region)
    except:
       raise(ValueError("Designated wildcard cannot be found in region list."))
        
    return rawRegionList[idx] 

def combineOutputs(build, sampleList_xr, sampleList_ds, sampleList_markers, regions="", outformat="real"):


    window = "_combined_"

    inputList = []
    if outformat == "real":
        for sample in sampleList_xr:
            sampledir = "results/XR/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + 
            "_plus_" + regions + window + "rpkm.txt")
            inputList.append(sampledir + sample + "_" + build + 
            "_minus_" + regions + window + "rpkm.txt")

        for sample in sampleList_ds:
            sampledir = "results/DS/" + sample + "/" 
            
            inputList.append(sampledir + sample + "_" + build + 
            "_plus_" + regions + window + "rpkm.txt")
            inputList.append(sampledir + sample + "_" + build + 
            "_minus_" + regions + window + "rpkm.txt")

    elif outformat == "intergenic":
        for sample in sampleList_xr:
            sampledir = "results/intergenic/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + 
            "_plus_" + regions + window + "rpkm.txt")
            inputList.append(sampledir + sample + "_" + build + 
            "_minus_" + regions + window + "rpkm.txt")

        for sample in sampleList_ds:
            sampledir = "results/intergenic/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + 
            "_plus_" + regions + window + "rpkm.txt")
            inputList.append(sampledir + sample + "_" + build + 
            "_minus_" + regions + window + "rpkm.txt")

    elif outformat == "sim":
        for sample in sampleList_xr:
            sampledir = "results/sim/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + 
            "_plus_" + regions + window + "rpkm.txt")
            inputList.append(sampledir + sample + "_" + build + 
            "_minus_" + regions + window + "rpkm.txt")

        for sample in sampleList_ds:
            sampledir = "results/sim/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + 
            "_plus_" + regions + window + "rpkm.txt")
            inputList.append(sampledir + sample + "_" + build + 
            "_minus_" + regions + window + "rpkm.txt")

    elif outformat == "sim_intergenic":
        for sample in sampleList_xr:
            sampledir = "results/intergenic_sim/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + 
            "_plus_" + regions + window + "rpkm.txt")
            inputList.append(sampledir + sample + "_" + build + 
            "_minus_" + regions + window + "rpkm.txt")

        for sample in sampleList_ds:
            sampledir = "results/intergenic_sim/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + 
            "_plus_" + regions + window + "rpkm.txt")
            inputList.append(sampledir + sample + "_" + build + 
            "_minus_" + regions + window + "rpkm.txt")  

    elif outformat == "markers_intergenic":
        for sample in sampleList_markers:
            sampledir = "results/intergenic_markers/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + regions + window + "rpkm.txt")

    elif outformat == "tss":
        for sample in sampleList_xr:
            sampledir = "results/XR/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + 
            "_sorted_tss_combined_rpkm.bed") 

        for sample in sampleList_ds:
            sampledir = "results/DS/" + sample + "/" 
            
            inputList.append(sampledir + sample + "_" + build + 
            "_sorted_dipyrimidines_tss_combined_rpkm.bed") 

    elif outformat == "tes":
        for sample in sampleList_xr:
            sampledir = "results/XR/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + 
            "_sorted_tes_combined_rpkm.bed") 

        for sample in sampleList_ds:
            sampledir = "results/DS/" + sample + "/" 
            
            inputList.append(sampledir + sample + "_" + build + 
            "_sorted_dipyrimidines_tes_combined_rpkm.bed") 

    #print(inputList)
    return inputList


def allInput(build="", sampleList=[], method="", regions=[]):

    inputList = []
    if method == "input":
    
        for sample in sampleList:
            sampledir = "results/input/" + sample + "/" 

            if isSingle(sample, sampleList, method):
                inputList.append(sampledir + sample + ".html")
            else:
                inputList.append(sampledir + sample + "_R1.html")
                inputList.append(sampledir + sample + "_R2.html")

            inputList.append(sampledir + sample + "_" + build + 
                "_length_distribution.png")
            inputList.append(sampledir + sample + "_" + build + 
                "_sorted_plus.bw")
            inputList.append(sampledir + sample + "_" + build + 
                "_sorted_minus.bw")

            inputList.append(sampledir + sample + "_" + 
                build + "_sorted_plus.bed") 
            inputList.append(sampledir + sample + "_" + 
                build + "_sorted_minus.bed") 
            inputList.append(sampledir + sample + "_" + 
                build + "_sim.bed") 

    if method == "okseq":
    
        for sample in sampleList:
            sampledir = "results/okseq/" + sample + "/" 

            if isSingle(sample, sampleList, method):
                inputList.append(sampledir + sample + ".html")
                inputList.append(sampledir + sample + "_se_" + build + 
                "_sorted.bam")
                inputList.append(sampledir + sample + "_se_" + build + 
                "_sorted.bam.bai")
            else:
                inputList.append(sampledir + sample + "_R1.html")
                inputList.append(sampledir + sample + "_R2.html")

    if method == "edu":

        for sample in sampleList:
            sampledir = "results/edu/" + sample + "/" 

            if isSingle(sample, sampleList, method):
                inputList.append(sampledir + sample + ".html")

            else:
                inputList.append(sampledir + sample + "_R1.html")
                inputList.append(sampledir + sample + "_R2.html")

            #if "inp" not in sample:
            #    inputList.append("results/edu/" + sample + "_" + build + 
            #    "_peaks.narrowPeak")
            #    inputList.append("results/edu/" + sample + "_" + build + 
            #    "_peaks.broadPeak")

            #inputList.append(sampledir + sample + "_" + build + "_sorted.bam")
            #inputList.append(sampledir + sample + "_" + build + 
            #    "_sorted_plus.bw")
            #inputList.append(sampledir + sample + "_" + build + 
            #    "_sorted_minus.bw")
            inputList.append("results/edu/merge_RT.txt")
            inputList.append("results/edu/merge_uv_RT.txt")
            inputList.append("results/edu/scatterplot_PearsonCorr_bigwigScores.png")
            inputList.append("results/edu/PearsonCorr_bigwigScores.tab")
            inputList.append("results/edu/heatmap_SpearmanCorr_readCounts.png")
            inputList.append("results/edu/SpearmanCorr_readCounts.tab")
            inputList.append("results/edu/PCA_readCounts.png")
            inputList.append("results/edu/PCA_readCounts_early.png")
            inputList.append("results/edu/PCA_readCounts_late.png")
            
    if method == "mutation":
    
        for sample in sampleList:
            sampledir = "results/mutation/" + sample + "/" 

            inputList.append(sampledir + sample + "_target_mut_plus.tsv") 
            inputList.append(sampledir + sample + "_target_mut_minus.tsv")

            for region in regions:
                inputList.append(sampledir + sample + "_target_mut_comb_" + 
                region + "_combined.txt")
                inputList.append(sampledir + sample + "_target_mut_comb_" + 
                region + "_intergenic_combined.txt")
                inputList.append(sampledir + sample + "_target_mut_comb_" + 
                region + "_org.txt")
                inputList.append(sampledir + sample + "_target_mut_comb_" + 
                region + "_intergenic_org.txt")
        
        for region in regions:

            inputList.append("results/mutation/" + region + "_counts.txt")

    if method == "ds":

        for sample in sampleList:
            sampledir = "results/DS/" + sample + "/" 
            simdir = "results/sim/" + sample + "/" 
            intdir = "results/intergenic/" + sample + "/" 
            intsimdir = "results/intergenic_sim/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + 
            "_sorted_dipyrimidines_tss_combined_rpkm.bed") 
            inputList.append(sampledir + sample + "_" + build + 
            "_sorted_dipyrimidines_tes_combined_rpkm.bed") 
            #inputList.append(sampledir + sample + "_" + build + 
            #"_sorted_dipyrimidines_TSNTS.bed")

            for region in regions:

                for mydir in [sampledir, simdir, intdir, intsimdir]:
                
                    inputList.append(mydir + sample + "_" + build + 
                    "_plus_" + region + "_combined_rpkm.txt")
                    inputList.append(mydir + sample + "_" + build + 
                    "_minus_" + region + "_combined_rpkm.txt")

    if method == "xr":

        for sample in sampleList:
            sampledir = "results/XR/" + sample + "/" 
            simdir = "results/sim/" + sample + "/" 
            intdir = "results/intergenic/" + sample + "/" 
            intsimdir = "results/intergenic_sim/" + sample + "/"

            inputList.append(sampledir + sample + "_" + build + 
            "_sorted_tss_combined_rpkm.bed") 
            inputList.append(sampledir + sample + "_" + build + 
            "_sorted_tes_combined_rpkm.bed") 
            #inputList.append(sampledir + sample + "_" + build + 
            #"_sorted_TSNTS.bed")

            for region in regions:

                for mydir in [sampledir, simdir, intdir, intsimdir]:

                    inputList.append(mydir + sample + "_" + build + 
                    "_plus_" + region + "_combined_rpkm.txt")
                    inputList.append(mydir + sample + "_" + build + 
                    "_minus_" + region + "_combined_rpkm.txt")

    if method == "markers":

        for sample in sampleList:

            sampledir = "results/intergenic_markers/" + sample + "/" 
            
            for region in regions:
        
                inputList.append(sampledir + sample + "_" + region + "_combined_rpkm.txt")

    if method == "report":

        inputList.append("results/plots/figure1B.pdf")
        inputList.append("results/plots/figure1C.pdf")
        inputList.append("results/plots/figure1D.pdf")
        inputList.append("results/plots/figure2.pdf")

        for region in regions:
                inputList.append("results/final/final_reports_" + build + 
                "_" + region + ".txt")
                inputList.append("results/final/final_reports_sim_" + build + 
                "_" + region + ".txt")
                inputList.append("results/final/final_reports_" + build + 
                "_" + region + "_intergenic.txt")
                inputList.append("results/final/final_reports_sim_" + build + 
                "_" + region + "_intergenic.txt")
                inputList.append("results/final/final_reports_markers_" + 
                region + "_intergenic.txt")

    #print(inputList)
    return inputList

################################################################################

#!/bin/env python 

import warnings
import os
import subprocess

################### Helper Functions ###########################################

def getSRR(sample, srrList, sampleList):

    try:
        idx = sampleList.index(sample)
    except:
       raise(ValueError("Designated wildcard cannot be found in sample list."))
        
    return srrList[idx]

def isSingle(sample, sampleList, srrEnabled, srrList, sample_dir):

    if srrEnabled:

        mySRR = getSRR(sample, srrList, sampleList)

        if mySRR == "NA":

            single = sample_dir + sample + ".fastq.gz"
            pairedR1 = sample_dir + sample + "_R1.fastq.gz"
            paired1 = sample_dir + sample + "_1.fastq.gz"
            
            if os.path.isfile(pairedR1) or os.path.isfile(paired1):
                return False
            elif os.path.isfile(single):
                return True
            else:
                raise(ValueError(paired1, single, "Sample not found..."))

        if ":" in mySRR:
            mySRR = mySRR.split(":")[0]

        shellCommand = 'fastq-dump -X 1 -Z --split-spot ' + mySRR + ' | wc -l'
        #print(shellCommand)
        p=subprocess.getoutput(shellCommand)
        #print(p)
        lineNum = int(p.split("\n")[2])
        #print(lineNum)

        if lineNum == 4:
            return True
        else:
            return False

    else:

        single = sample_dir + sample + ".fastq.gz"
        pairedR1 = sample_dir + sample + "_R1.fastq.gz"
        paired1 = sample_dir + sample + "_1.fastq.gz"
        
        if os.path.isfile(pairedR1) or os.path.isfile(paired1):
            return False
        elif os.path.isfile(single):
            return True
        else:
            raise(ValueError(paired1, single, "Sample not found..."))

def input4filter(wildcards, sampleList, srrEnabled, srrList, method, dirRaw):

    if isSingle(wildcards.samples, sampleList, srrEnabled, srrList, dirRaw):
        return "results/" + method + "/{samples}/{samples}_{build}_se.bed"
    else:    
        return "results/" + method + "/{samples}/{samples}_{build}_pe.bed"

def input4chip(wildcards, sampleList, srrEnabled, srrList):

    if isSingle(wildcards.samples, sampleList, srrEnabled, srrList, "resources/samples/chipseq/"):
        return "results/chipseq/{samples}/{samples}_{build}_se.bed"
    else:    
        return "results/chipseq/{samples}/{samples}_{build}_pe.bed"

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
                if line.strip().split(",")[0] == sample_name:

                    #targetLine = line.strip().replace(",","\t")
                    targetLine = line.strip()
                    return targetLine
                
            if targetLine == "":
                raise(ValueError(sample_name + " not found in markers.csv file..."))

def getMarkerName(wildcards, idx=1):
    return wildcards.samples.split("_")[idx]


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

            inputList.append(sampledir + sample + "_" + build + "_" + regions + window + "rpkm.txt")

    elif outformat == "methyl_intergenic":
        for sample in sampleList_markers:
            sampledir = "results/intergenic_methyl/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + "_" + regions + window + "rpkm.txt")

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


def allInput(build="", sampleList=[], srrEnabled=False, srrList=[], method="", regions=[]):

    inputList = []
    if method == "okseq":
    
        inputList.append("results/plots/PCA_readCounts_okseq.png")

        for sample in sampleList:
            sampledir = "results/okseq/" + sample + "/" 

            if isSingle(sample, sampleList, srrEnabled, srrList, "resources/samples/okseq/"):
                inputList.append(sampledir + sample + ".html")
                inputList.append(sampledir + sample + "_se_" + build + 
                "_sorted.bam")
                inputList.append(sampledir + sample + "_se_" + build + 
                "_sorted.bam.bai")
            else:
                inputList.append(sampledir + sample + "_R1.html")
                inputList.append(sampledir + sample + "_R2.html")
     
            inputList.append(sampledir + sample + "_" + build + "_HMMsegments_IZ.bed")

    if method == "edu":

        inputList.append("results/regions/repdomains_mean0.5.bed")
        inputList.append("results/regions/repdomains_uv_mean0.5.bed")
        inputList.append("results/plots/PCA_readCounts.png")
        inputList.append("results/plots/PCA_readCounts_early.png")
        inputList.append("results/plots/PCA_readCounts_late.png")

        for sample in sampleList:
            sampledir = "results/edu/" + sample + "/" 

            if isSingle(sample, sampleList, srrEnabled, srrList, "resources/samples/edu/"):
                inputList.append(sampledir + sample + ".html")

            else:
                inputList.append(sampledir + sample + "_R1.html")
                inputList.append(sampledir + sample + "_R2.html")

            inputList.append(sampledir + sample + "_" + build + 
                "_sorted_plus.bw")
            inputList.append(sampledir + sample + "_" + build + 
                "_sorted_minus.bw")
            
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
        
        inputList.append("results/plots/PCA_readCounts_chipseq.png")

        for sample in sampleList:

            sampledir = "results/intergenic_markers/" + sample + "/" 
            
            for region in regions:
        
                inputList.append(sampledir + sample + "_" + build + "_" + region + "_combined_rpkm.txt")

    if method == "methyl":

        for sample in sampleList:

            sampledir = "results/intergenic_methyl/" + sample + "/" 
            
            for region in regions:
        
                inputList.append(sampledir + sample + "_" + build + "_" + region + "_combined_rpkm.txt")

    if method == "report":
    
        inputList.append("results/plots/figure_markers.pdf")
        inputList.append("results/plots/figure1.pdf")
        inputList.append("results/plots/figure2.pdf")
        inputList.append("results/plots/figure3.pdf")
        inputList.append("results/plots/figure4.pdf")
        inputList.append("results/plots/figure5.pdf")
        inputList.append("results/plots/figureS2.pdf")
        inputList.append("results/plots/figureS3.pdf")
        inputList.append("results/plots/figureS4.pdf")

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
                #inputList.append("results/final/final_reports_methyl_" + 
                #region + "_intergenic.txt")

    #print(inputList)
    return inputList

################################################################################

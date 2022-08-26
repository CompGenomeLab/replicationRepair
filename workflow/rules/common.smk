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

            single = f"{sample_dir}{sample}.fastq.gz"
            pairedR1 = f"{sample_dir}{sample}_R1.fastq.gz"
            paired1 = f"{sample_dir}{sample}_1.fastq.gz"
            
            if os.path.isfile(pairedR1) or os.path.isfile(paired1):
                return False
            elif os.path.isfile(single):
                return True
            else:
                raise(ValueError(f"{paired1}, {pairedR1}, or {single} not found..."))

        if ":" in mySRR:
            mySRR = mySRR.split(":")[0]

        shellCommand = f"fastq-dump -X 1 -Z --split-spot {mySRR} | wc -l"
        
        p=subprocess.getoutput(shellCommand)
        
        lineNum = int(p.split("\n")[2])

        if lineNum == 4:
            return True
        else:
            return False

    else:

        single = f"{sample_dir}{sample}.fastq.gz"
        pairedR1 = f"{sample_dir}{sample}_R1.fastq.gz"
        paired1 = f"{sample_dir}{sample}_1.fastq.gz"
        
        if os.path.isfile(pairedR1) or os.path.isfile(paired1):
            return False
        elif os.path.isfile(single):
            return True
        else:
            raise(ValueError(f"{paired1}, {pairedR1}, or {single} not found..."))

def input4filter(wildcards, sampleList, srrEnabled, srrList, method, dirRaw):

    if isSingle(wildcards.samples, sampleList, srrEnabled, srrList, dirRaw):
        return "results/" + method + "/{samples}/{samples}_hg19_se.bed"
    else:    
        return "results/" + method + "/{samples}/{samples}_hg19_pe.bed"

def lineNum(file):
    
    linenum = 0
    if os.path.exists(file):
        with open(file) as f:
            for line in f:
                linenum += 1

        if linenum == 0:
            warnings.warn(f"\n{file} file is empty!")

    return linenum

def mappedReads(*files):

    lineNumber = 0
    for file in files:
        lineNumber += lineNum(str(file))

    return lineNumber

def info(wildcards):

    sample_name = wildcards.samples.replace("-","_")

    with open("resources/samples/samples.csv") as f:

        targetLine = ""
        for line in f:

            if line.strip().split(",")[1] == sample_name:

                #targetLine = line.strip().replace(",","\t")
                targetLine = line.strip()
                return targetLine
            
        if targetLine == "":
            raise(ValueError(f"{sample_name} not found in samples.csv file..."))

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

def combineOutputs(sampleList_xr, sampleList_ds, regions="", outformat="real"):

    inputList = []
    if outformat == "real":

        for sample in sampleList_xr:
            sampledir = f"results/XR/{sample}/" 

            inputList.append(f"{sampledir}{sample}_hg19_plus_{regions}_combined_rpkm.txt")
            inputList.append(f"{sampledir}{sample}_hg19_minus_{regions}_combined_rpkm.txt")

        for sample in sampleList_ds:
            sampledir = f"results/DS/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}_hg19_plus_{regions}_combined_rpkm.txt")
            inputList.append(f"{sampledir}{sample}_hg19_minus_{regions}_combined_rpkm.txt")

    elif outformat == "intergenic":

        for sample in (sampleList_xr + sampleList_ds):
            sampledir = f"results/intergenic/{sample}/" 

            inputList.append(f"{sampledir}{sample}_hg19_plus_{regions}_combined_rpkm.txt")
            inputList.append(f"{sampledir}{sample}_hg19_minus_{regions}_combined_rpkm.txt")

    elif outformat == "sim":

        for sample in (sampleList_xr + sampleList_ds):
            sampledir = f"results/sim/{sample}/"

            inputList.append(f"{sampledir}{sample}_hg19_plus_{regions}_combined_rpkm.txt")
            inputList.append(f"{sampledir}{sample}_hg19_minus_{regions}_combined_rpkm.txt")

    elif outformat == "sim_intergenic":

        for sample in (sampleList_xr + sampleList_ds):
            sampledir = f"results/intergenic_sim/{sample}/" 

            inputList.append(f"{sampledir}{sample}_hg19_plus_{regions}_combined_rpkm.txt")
            inputList.append(f"{sampledir}{sample}_hg19_minus_{regions}_combined_rpkm.txt")

    elif outformat == "tss":
        for sample in sampleList_xr:
            sampledir = f"results/XR/{sample}/" 

            inputList.append(f"{sampledir}{sample}_hg19_sorted_tss_combined_rpkm.bed") 

        for sample in sampleList_ds:
            sampledir = f"results/DS/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}_hg19_sorted_dipyrimidines_tss_combined_rpkm.bed") 

    elif outformat == "tes":
        for sample in sampleList_xr:
            sampledir = f"results/XR/{sample}/" 

            inputList.append(f"{sampledir}{sample}_hg19_sorted_tes_combined_rpkm.bed") 

        for sample in sampleList_ds:
            sampledir = f"results/DS/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}_hg19_sorted_dipyrimidines_tes_combined_rpkm.bed") 

    elif outformat == "tss_sim":
        for sample in sampleList_xr:
            sampledir = f"results/sim/{sample}/" 

            inputList.append(f"{sampledir}{sample}_hg19_sorted_tss_combined_rpkm.bed") 

        for sample in sampleList_ds:
            sampledir = f"results/sim/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}_hg19_sorted_dipyrimidines_tss_combined_rpkm.bed") 

    elif outformat == "tes_sim":
        for sample in sampleList_xr:
            sampledir = f"results/sim/{sample}/" 

            inputList.append(f"{sampledir}{sample}_hg19_sorted_tes_combined_rpkm.bed") 

        for sample in sampleList_ds:
            sampledir = f"results/sim/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}_hg19_sorted_dipyrimidines_tes_combined_rpkm.bed") 

    elif outformat == "tss_4kb":
        for sample in sampleList_xr:
            sampledir = f"results/XR/{sample}/" 

            inputList.append(f"{sampledir}{sample}_hg19_sorted_tss_windows_2001_1_combined_rpkm.bed") 

        for sample in sampleList_ds:
            sampledir = f"results/DS/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}_hg19_sorted_dipyrimidines_tss_windows_2001_1_combined_rpkm.bed") 

    elif outformat == "tes_4kb":
        for sample in sampleList_xr:
            sampledir = f"results/XR/{sample}/" 

            inputList.append(f"{sampledir}{sample}_hg19_sorted_tes_windows_2001_1_combined_rpkm.bed") 

        for sample in sampleList_ds:
            sampledir = f"results/DS/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}_hg19_sorted_dipyrimidines_tes_windows_2001_1_combined_rpkm.bed") 

    elif outformat == "tss_4kb_sim":
        for sample in sampleList_xr:
            sampledir = f"results/sim/{sample}/" 

            inputList.append(f"{sampledir}{sample}_hg19_sorted_tss_windows_2001_1_combined_rpkm.bed") 

        for sample in sampleList_ds:
            sampledir = f"results/sim/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}_hg19_sorted_dipyrimidines_tss_windows_2001_1_combined_rpkm.bed") 

    elif outformat == "tes_4kb_sim":
        for sample in sampleList_xr:
            sampledir = f"results/sim/{sample}/" 

            inputList.append(f"{sampledir}{sample}_hg19_sorted_tes_windows_2001_1_combined_rpkm.bed") 

        for sample in sampleList_ds:
            sampledir = f"results/sim/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}_hg19_sorted_dipyrimidines_tes_windows_2001_1_combined_rpkm.bed") 

    return inputList


def allInput(sampleList=[], srrEnabled=False, srrList=[], method="", regions=[]):

    inputList = []
    if method == "okseq":
    
        inputList.append("results/plots/PCA_readCounts_okseq.png")

        for sample in sampleList:
            sampledir = f"results/okseq/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}.html")

    if method == "edu":

        inputList.append("results/regions/repdomains_mean0.5.bed")
        inputList.append("results/regions/repdomains_uv_mean0.5.bed")
        inputList.append("results/plots/PCA_readCounts.png")
        inputList.append("results/plots/PCA_readCounts_early.png")
        inputList.append("results/plots/PCA_readCounts_late.png")

        for sample in sampleList:
            sampledir = f"results/edu/{sample}/"

            if isSingle(sample, sampleList, srrEnabled, srrList, "resources/samples/edu/"):
                inputList.append(f"{sampledir}{sample}.html")

            else:
                inputList.append(f"{sampledir}{sample}_R1.html")
                inputList.append(f"{sampledir}{sample}_R2.html")

    if method == "report":

        for fig_num in ["1", "2", "3", "4B", "4C_4D", "4E", "4F", "4G", "5A", 
        "5B_5C_5D", "S2A", "S3B_S3C_S3D", "S4", "S5", "S6", "S7", "S8", "S9", 
        "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19",
        "S20", "S21", "S22", "S23", "S24", "S25"]:
    
            inputList.append(f"results/plots/figure{fig_num}.pdf")

        for df_num in ["1_B1", "1_B2", "1_C1", "1_C2", "1_C3", "1_C4", "1_D",
        "2_A1", "2_A2", "2_B", "3_A", "3_B", "4_B1", "4_B2", "4_B3", "4_B4", 
        "4_C1", "4_C2", "4_D1", "4_D2", "4_E", "4_F", "4_G", "5_B", "5_C", 
        "5_D1", "5_D2", "S3_B", "S3_C1", "S3_C2", "S3_D", "S4", "S5_A", "S5_B", 
        "S6_A", "S6_B", "S7_A", "S7_B", "S8", "S9_A1", "S9_A2", "S9_B1", 
        "S9_B2", "S10_A1", "S10_A2", "S10_B1", "S10_B2", "S11_A1", "S11_A2", 
        "S11_B1", "S11_B2", "S12_A1", "S12_A2", "S12_B1", "S12_B2", "S13_A1", 
        "S13_A2", "S13_B1", "S13_B2", "S14", "S15_A1", "S15_A2", "S15_B1", 
        "S15_B2", "S16_B1", "S16_B2", "S16_C1", "S16_C2", "S17", "S18", "S19", 
        "S20", "S21", "S22", "S23", "S24", "S25"]:

            inputList.append(f"results/plot_dataframe/figure{df_num}.csv")

    return inputList

################################################################################

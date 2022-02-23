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

    warnMessage = (f"\n{file} file is either empty or does not exists!\n" + 
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
                raise(ValueError(f"{sample_name} not found in samples.csv file..."))
    
    elif method == "marker":

        with open("resources/samples/markers.csv") as f:

            targetLine = ""
            for line in f:
                if line.strip().split(",")[0] == sample_name:

                    #targetLine = line.strip().replace(",","\t")
                    targetLine = line.strip()
                    return targetLine
                
            if targetLine == "":
                raise(ValueError(f"{sample_name} not found in markers.csv file..."))

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

    inputList = []
    if outformat == "real":

        for sample in sampleList_xr:
            sampledir = f"results/XR/{sample}/" 

            inputList.append(f"{sampledir}{sample}_{build}_plus_{regions}_combined_rpkm.txt")
            inputList.append(f"{sampledir}{sample}_{build}_minus_{regions}_combined_rpkm.txt")

        for sample in sampleList_ds:
            sampledir = f"results/DS/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}_{build}_plus_{regions}_combined_rpkm.txt")
            inputList.append(f"{sampledir}{sample}_{build}_minus_{regions}_combined_rpkm.txt")

    elif outformat == "intergenic":

        for sample in (sampleList_xr + sampleList_ds):
            sampledir = f"results/intergenic/{sample}/" 

            inputList.append(f"{sampledir}{sample}_{build}_plus_{regions}_combined_rpkm.txt")
            inputList.append(f"{sampledir}{sample}_{build}_minus_{regions}_combined_rpkm.txt")

    elif outformat == "sim":

        for sample in (sampleList_xr + sampleList_ds):
            sampledir = f"results/sim/{sample}/"

            inputList.append(f"{sampledir}{sample}_{build}_plus_{regions}_combined_rpkm.txt")
            inputList.append(f"{sampledir}{sample}_{build}_minus_{regions}_combined_rpkm.txt")

    elif outformat == "sim_intergenic":

        for sample in (sampleList_xr + sampleList_ds):
            sampledir = f"results/intergenic_sim/{sample}/" 

            inputList.append(f"{sampledir}{sample}_{build}_plus_{regions}_combined_rpkm.txt")
            inputList.append(f"{sampledir}{sample}_{build}_minus_{regions}_combined_rpkm.txt")

    elif outformat == "markers_intergenic":

        for sample in sampleList_markers:
            sampledir = f"results/intergenic_markers/{sample}/" 

            inputList.append(f"{sampledir}{sample}_{build}_{regions}_combined_rpkm.txt")

    elif outformat == "methyl_intergenic":

        for sample in sampleList_markers:
            sampledir = f"results/intergenic_methyl/{sample}/"  

            inputList.append(f"{sampledir}{sample}_{build}_plus_{regions}_combined_rpkm.txt")
            inputList.append(f"{sampledir}{sample}_{build}_minus_{regions}_combined_rpkm.txt")

    elif outformat == "tss":
        for sample in sampleList_xr:
            sampledir = f"results/XR/{sample}/" 

            inputList.append(f"{sampledir}{sample}_{build}_sorted_tss_combined_rpkm.bed") 

        for sample in sampleList_ds:
            sampledir = f"results/DS/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}_{build}_sorted_dipyrimidines_tss_combined_rpkm.bed") 

    elif outformat == "tes":
        for sample in sampleList_xr:
            sampledir = f"results/XR/{sample}/" 

            inputList.append(f"{sampledir}{sample}_{build}_sorted_tes_combined_rpkm.bed") 

        for sample in sampleList_ds:
            sampledir = f"results/DS/{sample}/" 
            
            inputList.append(f"{sampledir}{sample}_{build}_sorted_dipyrimidines_tes_combined_rpkm.bed") 

    #print(inputList)
    return inputList


def allInput(build="", sampleList=[], srrEnabled=False, srrList=[], method="", regions=[]):

    inputList = []
    if method == "okseq":
    
        inputList.append("results/plots/PCA_readCounts_okseq.png")

        for sample in sampleList:
            sampledir = f"results/okseq/{sample}/" 

            if isSingle(sample, sampleList, srrEnabled, srrList, "resources/samples/okseq/"):
                inputList.append(f"{sampledir}{sample}.html")
                inputList.append(f"{sampledir}{sample}_se_{build}_sorted.bam")
                inputList.append(f"{sampledir}{sample}_se_{build}_sorted.bam.bai")
            else:
                inputList.append(f"{sampledir}{sample}_R1.html")
                inputList.append(f"{sampledir}{sample}_R2.html")
     
            inputList.append(f"{sampledir}{sample}_{build}_HMMsegments_IZ.bed")

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

            inputList.append(f"{sampledir}{sample}_{build}_sorted_plus.bw")
            inputList.append(f"{sampledir}{sample}_{build}_sorted_minus.bw")
            
    if method == "mutation":
    
        for sample in sampleList:
            sampledir = f"results/mutation/{sample}/" 
            sampledir_int = f"results/intergenic_mutation/{sample}/"

            inputList.append(f"{sampledir}{sample}_target_mut_plus.tsv") 
            inputList.append(f"{sampledir}{sample}_target_mut_minus.tsv")
            inputList.append(f"{sampledir}{sample}_intergenic_target_mut_plus.tsv") 
            inputList.append(f"{sampledir}{sample}_intergenic_target_mut_minus.tsv")

            for region in regions:
                inputList.append(f"{sampledir}{sample}_target_mut_{region}_combined_rpkm.txt")
                inputList.append(f"{sampledir_int}{sample}_target_mut_{region}_combined_rpkm.txt")
        
        for region in regions:

            inputList.append(f"results/regions/{region}_counts.txt")

    if method == "ds":

        for sample in sampleList:
            sampledir = f"results/DS/{sample}/" 
            simdir = f"results/sim/{sample}/" 
            intdir = f"results/intergenic/{sample}/" 
            intsimdir = f"results/intergenic_sim/{sample}/" 

            inputList.append(f"{sampledir}{sample}_{build}_sorted_dipyrimidines_tss_combined_rpkm.bed") 
            inputList.append(f"{sampledir}{sample}_{build}_sorted_dipyrimidines_tes_combined_rpkm.bed") 
            #inputList.append(f"{sampledir}{sample}_{build}_sorted_dipyrimidines_TSNTS.bed")

            for region in regions:

                for mydir in [sampledir, simdir, intdir, intsimdir]:
                
                    inputList.append(f"{mydir}{sample}_{build}_plus_{region}_combined_rpkm.txt")
                    inputList.append(f"{mydir}{sample}_{build}_minus_{region}_combined_rpkm.txt")

    if method == "xr":

        for sample in sampleList:
            sampledir = f"results/XR/{sample}/" 
            simdir = f"results/sim/{sample}/" 
            intdir = f"results/intergenic/{sample}/" 
            intsimdir = f"results/intergenic_sim/{sample}/"

            inputList.append(f"{sampledir}{sample}_{build}_sorted_tss_combined_rpkm.bed") 
            inputList.append(f"{sampledir}{sample}_{build}_sorted_tes_combined_rpkm.bed") 
            #inputList.append(f"{sampledir}{sample}_{build}_sorted_TSNTS.bed")

            for region in regions:

                for mydir in [sampledir, simdir, intdir, intsimdir]:

                    inputList.append(f"{mydir}{sample}_{build}_plus_{region}_combined_rpkm.txt")
                    inputList.append(f"{mydir}{sample}_{build}_minus_{region}_combined_rpkm.txt")

    if method == "markers":
        
        inputList.append("results/plots/PCA_readCounts_chipseq.png")

        for sample in sampleList:

            sampledir = f"results/intergenic_markers/{sample}/" 
            
            for region in regions:
        
                inputList.append(f"{sampledir}{sample}_{build}_{region}_combined_rpkm.txt")

    if method == "methyl":

        for sample in sampleList:

            sampledir = "results/intergenic_methyl/{sample}/" 
            
            for region in regions:
        
                inputList.append(f"{sampledir}{sample}_{build}_minus_{region}_combined_rpkm.txt")
                inputList.append(f"{sampledir}{sample}_{build}_plus_{region}_combined_rpkm.txt")

    if method == "report":
    
        inputList.append("results/plots/figure1.pdf")
        inputList.append("results/plots/figure2.pdf")
        inputList.append("results/plots/figure3.pdf")
        inputList.append("results/plots/figure4B.pdf")
        inputList.append("results/plots/figure4C_4D.pdf")
        inputList.append("results/plots/figure4E.pdf")
        inputList.append("results/plots/figure5A.pdf")
        inputList.append("results/plots/figure5B_5C_5D.pdf")
        inputList.append("results/plots/figureS2A.png")
        inputList.append("results/plots/figureS3B_S3C_S3D.pdf")
        inputList.append("results/plots/figureS4.pdf")
        inputList.append("results/plots/figureS5.pdf")
        inputList.append("results/plots/figureS6.pdf")
        inputList.append("results/plots/figureS7.pdf")
        inputList.append("results/plots/figureS8.pdf")
        inputList.append("results/plots/figureS9.pdf")
        inputList.append("results/plots/figureS10.pdf")
        inputList.append("results/plots/figureS11.pdf")
        inputList.append("results/plots/figureS12B_S12C.pdf")
        inputList.append("results/plots/figureS13.pdf")
        inputList.append("results/plots/figureS14.pdf")
        inputList.append("results/plots/figureS15.pdf")
        inputList.append("results/plots/figureS16.pdf")
        inputList.append("results/plots/figureS17.pdf")


        #inputList.append("results/plots/figure_markers.pdf")
        #inputList.append("results/plots/figure_methyl.pdf")
        inputList.append("results/plots/figure_normDSXR.pdf")

        for region in regions:
                inputList.append(f"results/final/final_reports_{build}_{region}.txt")
                inputList.append(f"results/final/final_reports_sim_{build}_{region}.txt")
                inputList.append(f"results/final/final_reports_{build}_{region}_intergenic.txt")
                inputList.append(f"results/final/final_reports_sim_{build}_{region}_intergenic.txt")
                #inputList.append(f"results/final/final_reports_markers_{region}_intergenic.txt")
                #inputList.append(f"results/final/final_reports_methyl_{region}_intergenic.txt")

    return inputList

################################################################################

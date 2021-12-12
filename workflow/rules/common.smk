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

def getInput(sample, inputExist, inputList, inputIdx, sampleList, build, region=""):

    if inputExist:
        inpDict={}
        for inp_idx in range(len(inputIdx)):
            idx_split = inputIdx[inp_idx].strip().split(",")
            indexList=[]
            for sample_idx in idx_split:
                sample_idx = sample_idx.strip() 
                if "-" in sample_idx:
                    for range_idx in range(int(sample_idx.split("-")[0]), int(sample_idx.split("-")[1])+1):
                        indexList.append(int(range_idx)) 
                else:
                    indexList.append(int(sample_idx)) 
                
            for sample_idx in indexList:
                if inputList[inp_idx] not in inpDict:
                    inpDict[inputList[inp_idx]] = [sampleList[sample_idx]]
                else:
                    inpDict[inputList[inp_idx]].append(sampleList[sample_idx])
        for k,v in inpDict.items():
        
            if sample in v and region=="IZ":
            
                return "resources/samples/input/" + k + "_" + build + "_IZ.fasta"
            
            if sample in v:
            
                return "resources/samples/input/" + k + "_" + build + ".fasta"
                
    else:
        return "resources/ref_genomes/" + build + "/genome_" + build + ".ron" 


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

def combineOutputs(build, sampleList_xr, sampleList_ds, sampleList_markers, regions="", outformat="real", kmer=""):


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

    elif outformat == "sim_kmer":

        for sample in sampleList_xr:
            sampledir = "results/sim/" + sample + "/" 

            inputList.append(sampledir + sample + "_xr_" + kmer + "_" + build + 
            "_plus_" + regions + window + "rpkm.txt")
            inputList.append(sampledir + sample + "_xr_" + kmer + "_" + build + 
            "_minus_" + regions + window + "rpkm.txt")

        for sample in sampleList_ds:
            sampledir = "results/sim/" + sample + "/" 

            inputList.append(sampledir + sample + "_ds_" + kmer + "_" + build + 
            "_plus_" + regions + window + "rpkm.txt")
            inputList.append(sampledir + sample + "_ds_" + kmer + "_" + build + 
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

            inputList.append(sampledir + sample + "_" + build + "_minus_" + regions + window + "rpkm.txt")
            inputList.append(sampledir + sample + "_" + build + "_plus_" + regions + window + "rpkm.txt")

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
            sampledir_int = "results/intergenic_mutation/" + sample + "/" 

            inputList.append(sampledir + sample + "_target_mut_plus.tsv") 
            inputList.append(sampledir + sample + "_target_mut_minus.tsv")
            inputList.append(sampledir + sample + "_intergenic_target_mut_plus.tsv") 
            inputList.append(sampledir + sample + "_intergenic_target_mut_minus.tsv")

            for region in regions:
                inputList.append(sampledir + sample + "_target_mut_" + 
                region + "_combined_rpkm.txt")
                inputList.append(sampledir_int + sample + "_target_mut_" + 
                region + "_combined_rpkm.txt")
        
        for region in regions:

            inputList.append("results/regions/" + region + "_counts.txt")

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
        
                inputList.append(sampledir + sample + "_" + build + "_minus_" + region + "_combined_rpkm.txt")
                inputList.append(sampledir + sample + "_" + build + "_plus_" + region + "_combined_rpkm.txt")

    if method == "report":
    
        inputList.append("results/regions/genome_T_counts.txt")
        inputList.append("results/plots/figure1.pdf")
        inputList.append("results/plots/figure2.pdf")
        inputList.append("results/plots/figure3.pdf")
        inputList.append("results/plots/figure4.pdf")
        inputList.append("results/plots/figure5.pdf")
        inputList.append("results/plots/figure6A.pdf")
        inputList.append("results/plots/figure6.pdf")
        inputList.append("results/plots/figureS2.pdf")
        inputList.append("results/plots/figureS3.pdf")
        inputList.append("results/plots/figureS4.pdf")
        inputList.append("results/plots/figureS5.pdf")
        inputList.append("results/plots/figureS5_intergenic.pdf")
        inputList.append("results/plots/figureS5_64.pdf")
        inputList.append("results/plots/figureS5_intergenic_64.pdf")  
        inputList.append("results/plots/figureS5_repdomains.pdf")
        inputList.append("results/plots/figureS5_repdomains_intergenic.pdf")
        inputList.append("results/plots/figureS5_repdomains_64.pdf")
        inputList.append("results/plots/figureS5_repdomains_intergenic_64.pdf")  
        #inputList.append("results/plots/figureS5_repdomains_2.pdf")
        #inputList.append("results/plots/figureS5_repdomains_3.pdf")
        #inputList.append("results/plots/figureS5_repdomains_4.pdf")
        #inputList.append("results/plots/figureS5_repdomains_5.pdf")
        #inputList.append("results/plots/figureS5_repdomains_rmTTTT.pdf")
        #inputList.append("results/plots/figureS5_repdomains_IZ_rmTTTT.pdf")
        #inputList.append("results/plots/figureS5_repdomains_IZ.pdf")
        inputList.append("results/plots/figureS6.pdf")
        inputList.append("results/plots/figureS7.pdf")
        inputList.append("results/plots/figureS8.pdf")
        inputList.append("results/plots/figureS9.pdf")
        inputList.append("results/plots/figure_markers.pdf")
        inputList.append("results/plots/figure_methyl.pdf")
        inputList.append("results/plots/mer1.pdf")
        inputList.append("results/plots/mer2.pdf")
        inputList.append("results/plots/mer3.pdf")
        inputList.append("results/plots/mer4.pdf")
        inputList.append("results/plots/mer5.pdf")
        inputList.append("results/plots/mer.pdf")
        inputList.append("results/table/seq_asymmetry.csv")

        for region in regions:
                inputList.append("results/final/final_reports_" + build + 
                "_" + region + ".txt")
                inputList.append("results/final/final_reports_sim_" + build + 
                "_" + region + ".txt")
                #inputList.append("results/final/final_reports_sim_2_" + build + 
                #"_" + region + ".txt")
                #inputList.append("results/final/final_reports_sim_3_" + build + 
                #"_" + region + ".txt")
                #inputList.append("results/final/final_reports_sim_4_" + build + 
                #"_" + region + ".txt")
                #inputList.append("results/final/final_reports_sim_5_" + build + 
                #"_" + region + ".txt")
                inputList.append("results/final/final_reports_sim_IZ_" + build + 
                "_" + region + ".txt")
                inputList.append("results/final/final_reports_sim_IZ_rmTTTT_" + build + 
                "_" + region + ".txt")
                #inputList.append("results/final/final_reports_sim_rmTTTT_" + build + 
                #"_" + region + ".txt")
                inputList.append("results/final/final_reports_" + build + 
                "_" + region + "_intergenic.txt")
                inputList.append("results/final/final_reports_sim_" + build + 
                "_" + region + "_intergenic.txt")
                inputList.append("results/final/final_reports_markers_" + 
                region + "_intergenic.txt")
                inputList.append("results/final/final_reports_methyl_" + 
                region + "_intergenic.txt")

    #print(inputList)
    return inputList

################################################################################

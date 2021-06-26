
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
        raise(ValueError("Sample not found..."))


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

def info(wildcards):
    sample_name = wildcards.samples.replace("-","_")

    with open("resources/samples/samples.csv") as f:

        targetLine = ""
        for line in f:

            if line.split(",")[1] == sample_name:

                targetLine = line.replace(",","\t")
                return targetLine
            
        if targetLine == "":
            raise(ValueError(sample_name + " not found in samples.csv file..."))

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

def allInput(build, sampleList, method, regions=[]):

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

    if method == "mutation":
    
        for sample in sampleList:
            sampledir = "results/mutation/" + sample + "/" 

            inputList.append(sampledir + sample + "_target_mut_plus.tsv") 
            inputList.append(sampledir + sample + "_target_mut_minus.tsv")

    if method == "ds":

        for sample in sampleList:
            sampledir = "results/DS/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + 
            "_sorted_dipyrimidines_tss_combined_rpkm.bed") 
            inputList.append(sampledir + sample + "_" + build + 
            "_sorted_dipyrimidines_TSNTS.bed")

            for region in regions:
                inputList.append(sampledir + sample + "_" + build + 
                "_sorted_ds_dipyrimidines_plus_" + region + "_combined_rpkm.txt")
                inputList.append(sampledir + sample + "_" + build + 
                "_sorted_ds_dipyrimidines_minus_" + region + "_combined_rpkm.txt")

    if method == "xr":

        for sample in sampleList:
            sampledir = "results/XR/" + sample + "/" 

            inputList.append(sampledir + sample + "_" + build + 
            "_sorted_tss_combined_rpkm.bed") 
            inputList.append(sampledir + sample + "_" + build + 
            "_sorted_TSNTS.bed")

            for region in regions:
                inputList.append(sampledir + sample + "_" + build + 
                "_sorted_xr_plus_" + region + "_combined_rpkm.txt")
                inputList.append(sampledir + sample + "_" + build + 
                "_sorted_xr_minus_" + region + "_combined_rpkm.txt")

    #print(inputList)
    return inputList

################################################################################

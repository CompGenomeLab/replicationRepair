
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

def input4filter(wildcards, sampleList, method):

    if isSingle(wildcards.samples, sampleList, method):
        return "results/" + method + "/{samples}/{samples}_{build}_se.bed"
    else:    
        return "results/" + method + "/{samples}/{samples}_{build}_pe.bed"

def input4nucTable(method):

    if method == "XR":
        return "results/{samples}/{samples}_{build}_lengthMode.fa"
    elif method == "DS":    
        return "results/{samples}/{samples}_{build}_sorted_10.fa"

def getMotif(wildcards):
    
    if "oxaliplatin" in wildcards.samples or "cisplatin" in wildcards.samples: 
        return "'.{4}(g|G){2}.{4}'"
    
    elif "64" in wildcards.samples or "CPD" in wildcards.samples:
        return "'.{4}(c|t|C|T){2}.{4}'"

def getDinuc(wildcards):
    
    if "Oxaliplatin" in wildcards.samples or "Cisplatin" in wildcards.samples: 
        return "'GG'"
    
    elif "64" in wildcards.samples or "CPD" in wildcards.samples or "R190" in wildcards.samples:
        return "'CC','CT','TC','TT'"

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

        for line in f:

            if line.split(",")[1] == sample_name:

                return line.replace(",","\t")
            
            else:
                raise(ValueError(sample_name + " not found in samples.csv file..."))

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
                "_sorted_ds_dipyrimidines_plus_" + region + ".txt")

    #print(inputList)
    return inputList

################################################################################

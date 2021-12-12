#!/usr/bin/env python

import os
import sys
import argparse


parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='fasta to bed')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input fasta (stranded)')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output bed')
args = parser.parse_args()

filein = args.i
out = args.o


def fa2list(ifile):

    bed_list=[]
    data_line = ifile.readlines()
    for line in data_line:

        if line[0] == ">":
            strand = line.strip()[-2]
            line = line.replace(":","-").replace("(","-").replace(">","-")
            lineList = line.split("-")  
            bed = "\t".join([lineList[1],lineList[2],lineList[3],".",".",strand])
            bed_list.append(bed)
    
    return bed_list

def list2bed(bed_list, out):

    for bed in bed_list:
        out.write(bed + "\n")

bedList = fa2list(filein)
list2bed(bedList, out)

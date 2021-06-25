#!/usr/bin/env python

import os
import sys
import argparse
import fileinput
import re
from collections import Counter

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='calculate RPKM values')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-c', nargs="+", required= True, default=sys.stdin, help='line number of the counts')
parser.add_argument('-mr', nargs="+", default=sys.stdin, help='line number of the mapped reads')
parser.add_argument('-chse', nargs="+", default=sys.stdin, help='line number of the chromosome start and end position respectively')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
args = parser.parse_args()

countline = (int((str(args.c).strip().split("'"))[1]) -1) 
#print(countline)
   
filein = args.i
out = args.o
separator = "\t"

for line in filein:

    #print(line)
    ll = line.strip().split(separator)
    chr_endline = int((str(args.chse).strip().split("'"))[3]) -1
    chr_startline = int((str(args.chse).strip().split("'"))[1]) -1        
    domain_length = float(ll[chr_endline]) - float(ll[chr_startline])
    #print(domain_length) 
    
    if type(args.mr) is list:
        mappedreadsline = (int((str(args.mr).strip().split("'"))[1]) -1)
        mappedReads = float(ll[mappedreadsline])
        #print(mappedreadsline)
    else: mappedReads = 1000000     

    RPKM = (1000000 * float(ll[countline]) / mappedReads) * 1000 / domain_length
    RPKM = [str(RPKM)]
    #print(RPKM)
    #print(float(ll[mappedreadsline]))
    newList = ll + RPKM
    #print(newList)
    out.write(separator.join(newList) + "\n")
      
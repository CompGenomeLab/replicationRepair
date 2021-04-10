#!/usr/bin/env python

import os
import sys
import argparse
import fileinput
from collections import Counter

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', default=sys.stdin, help='input')
parser.add_argument('-g', default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
args = parser.parse_args()

filein = args.i
genome = args.g
out = args.o
separator = "\t"

num_linesg = 0
with open(genome, 'r') as fg:
    for line in fg:
        num_linesg += 1
#print(num_linesg)

num_linesin = 0
with open(filein, 'r') as fin:
    for line in fin:
            num_linesin += 1
#print(num_linesin)

with open(genome, 'r') as fg:
        t=[x.strip().split(separator) for x in fg.readlines()]
        #print(t)
        max_int =[int(x[2]) for x in t]
        #print(max_int)
        min_int = [int(x[1]) for x in t]
        #print(min_int)
        ref_chr = [x[0] for x in t]
        #print(ref_chr)

with open(filein, 'r') as fin:        
    c=[y.strip().split(separator) for y in fin.readlines()]
    line = [y for y in c]
    #print(line)
    chr = [y[0] for y in c]
    #print(chr)
    mn = [int(y[1]) for y in c]   
    mx = [int(y[2]) for y in c]    
    #print(mn)
    #print(mx)

matchedLine = []
for i in range(0,num_linesg-1):
        for h in range(0,num_linesin-1):
                if ref_chr[i] == chr[h]:
                        if min_int[i] < mn[h] < mx[h] < max_int[i]: 
                            matchedLine = line[h] 
                            out.write(separator.join(matchedLine) + "\n")
                        
        
     




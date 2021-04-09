#!/usr/bin/env python

import os
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='count motifs for each header')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-m', '--motif', nargs="+", required= True, default=sys.stdin, help='provide motifs in a list')
parser.add_argument('--agg', action='store_true', help='aggregate positions by name')
args = parser.parse_args()

filein = args.i
out = args.o
motif = args.motif
separator = "\t"

def countMotif(seqLine, motif):

    return seqLine.seq.lower().count_overlap(motif.lower())

headerList = ["id"]
for m in motif:
    headerList.append(m)

out.write(separator.join(headerList) + "\n")

if args.agg:

    defaultDict = {'num': 0}
    for m in motif:
        defaultDict[m] = 0
        
    headerDict={}
    for r in SeqIO.parse(filein, 'fasta'):
        
        if r.name not in headerDict:
            headerDict[r.name] = defaultDict.copy()
            
        headerDict[r.name]['num'] += 1 
        for m in motif:
            headerDict[r.name][m] += countMotif(r, m)         

    for k,v in headerDict.items():
        
        outLine = [k]
        
        nestedDict = v.copy()
        del nestedDict['num']
        for m in nestedDict:
            outLine.append(str(nestedDict[m]/v['num']))

        out.write(separator.join(outLine))
        out.write("\n")
    
else:
    for r in SeqIO.parse(filein, 'fasta'):

        outLine = [r.name]
        for m in motif:
            outLine.append(str(countMotif(r, m)))

        out.write(separator.join(outLine))
        out.write("\n")

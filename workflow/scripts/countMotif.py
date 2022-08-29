#!/usr/bin/env python

import os
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='count motifs for each position')
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

if args.agg:

    headerList = [m for m in ["id"] + motif]

    out.write(separator.join(headerList) + "\n")

    defaultDict = {'num': 0}
    for m in motif:
        defaultDict[m] = 0
        
    headerDict={}
    for r in SeqIO.parse(filein, 'fasta'):
        
        if r.name.strip().split(":")[0] not in headerDict:
            headerDict[r.name.strip().split(":")[0]] = defaultDict.copy()
            
        headerDict[r.name.strip().split(":")[0]]['num'] += 1 
        for m in motif:
            headerDict[r.name.strip().split(":")[0]][m] += countMotif(r, m)         

    for k,v in headerDict.items():
        
        outLine = [k]
        
        nestedDict = v.copy()
        del nestedDict['num']
        for m in nestedDict:
            outLine.append(str(nestedDict[m]/v['num']))

        out.write(separator.join(outLine))
        out.write("\n")
    
else:
    headerList = [m for m in ["chr","start","end"] + motif]

    out.write(separator.join(headerList) + "\n")

    for r in SeqIO.parse(filein, 'fasta'):

        r.name = r.name.replace("-",":")
        outLine = r.name.strip().split(":")
        for m in motif:

            length = int(outLine[2]) - int(outLine[1])
            outLine.append(str(countMotif(r, m)/(length-len(m)+1)*100))

        out.write(separator.join(outLine))
        out.write("\n")

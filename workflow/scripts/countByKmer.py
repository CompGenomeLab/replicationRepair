#!/usr/bin/env python

import os
import sys
import argparse
from Bio import SeqIO
from more_itertools import distinct_permutations

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='count motifs for each header')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-m', '--kmer', nargs="+", required= True, default=sys.stdin, help='number of kmer')
parser.add_argument('--agg', action='store_true', help='aggregate positions by name')
args = parser.parse_args()

filein = args.i
out = args.o
kmer = int(args.kmer[0])
separator = "\t"

def kmerList(kmer):
    
    letters=""
    seqList=[]
    for i in range(1, kmer+1):
        letters+="A"
        letters+="T"
        letters+="C"
        letters+="G"

        for p in distinct_permutations(letters, i):
            seqList.append(''.join(p))

    return seqList

def countMotif(seqLine, motif):

    return seqLine.seq.lower().count_overlap(motif.lower())

motif = kmerList(kmer)
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
    for r in SeqIO.parse(filein, 'fasta'):

        outLine = [r.name.strip().split(":")[0]]
        for m in motif:
            outLine.append(str(countMotif(r, m)))

        out.write(separator.join(outLine))
        out.write("\n")

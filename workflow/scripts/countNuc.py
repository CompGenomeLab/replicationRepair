import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='count motif throughout genome')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-m', nargs='+', default=sys.stdin, help='motif to be calculated')
args = parser.parse_args()

filein = args.i
motif = args.m[0]

def countNuc(filein):

    data_line = filein.readlines()
    countChrom = 0
    genome = ""

    for line in data_line:

        if line[0] != ">":
            genome += line.strip()

        else:
            countChrom += len(motif)-1

    genome = genome.replace("N","")
  
    countTarget = genome.upper().count(motif)

    return countTarget/(len(genome)-countChrom)*100

print(countNuc(filein))
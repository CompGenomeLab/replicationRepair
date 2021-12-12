import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='count T throughout genome')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
args = parser.parse_args()

filein = args.i

def countT(filein):

    data_line = filein.readlines()
    countT = 0
    genome = ""

    for line in data_line:

        if line[0] != ">":
            genome += line.strip()
            
    countT = genome.upper().count("T")
    
    return countT/len(genome)*100

print(countT(filein))
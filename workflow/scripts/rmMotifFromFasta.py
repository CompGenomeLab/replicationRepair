import argparse
import sys

parser = argparse.ArgumentParser(description='converts fasta to bed by choosing a motif of interest')
parser.add_argument('-i', required= False, help='input')
parser.add_argument('-o', required= False, help='output')
parser.add_argument('-m', required= False, help='method that produced the reads (XR/DS)')
parser.add_argument('-r', required= True, help='regex motif that is expected')

args = parser.parse_args()

fastain = args.i
output = args.o
pattern = args.r
method = args.m

def fastaHeader2bedLine(header):
    chr = header.strip().replace('>','').split(':')[0]
    start = int(header.strip().split(':')[1].split('(')[0].split('-')[0])
    end = int(header.strip().split(':')[1].split('(')[0].split('-')[1])
    strand = header.strip().split(':')[1].split('(')[1][0]
    return chr + '\t' + str(start) + '\t' + str(end) + '\t.\t.\t' + strand

out = open(output, 'w')

with open (fastain) as fin:

    if method == "DS":

        for line in fin:

            if ">" in line:
                header = line

            else:
                if pattern not in line.strip()[3:7]:
                    out.write(fastaHeader2bedLine(header) + "\n")
    
    elif method == "XR":

        for line in fin:

            if ">" in line:
                header = line

            else:
                if pattern not in line.strip()[-9:-5]:
                    out.write(fastaHeader2bedLine(header) + "\n")     

out.close()
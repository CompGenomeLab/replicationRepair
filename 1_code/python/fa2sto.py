from Bio import SeqIO
import sys
import argparse


parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='collapse windows')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
args = parser.parse_args()

filein = args.i
out = args.o

records = SeqIO.parse(filein, "fasta")
count = SeqIO.write(records, out, "stockholm")
print("Converted %i records" % count)
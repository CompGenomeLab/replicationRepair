#!/usr/bin/env python

import os
import sys
import argparse


parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='collapse windows')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-chr', default=sys.stdout, help='keep chromosomes separate. Usage: "-chr T".')
parser.add_argument('-score', default=sys.stdout, help='keep scores separate. Use "-score T".')
parser.add_argument('-strand', default=sys.stdout, help='keep strands separate. Use "-strand T".')
args = parser.parse_args()

filein = args.i
out = args.o


def bedIntersected_dictionary(file=filein,chr=".",score="0",strand="."):

    bed_dict={ }
    data_line = file.readlines()
    for x in data_line:
        data_list = x.split()

        if args.chr == "T":
            chr = data_list[0]

        if args.score == "T": 
            score = str(data_list[4])   
        
        if args.strand == "T": 
            strand = data_list[5]

        length = int(data_list[2]) - int(data_list[1])
        name = data_list[3]            
        intersect = float(data_list[6])

        bed_key = "%s_cut_%s_cut_%s_cut_%s" %(name,chr,score,strand) 

        if bed_key in bed_dict:
            bed_dict[bed_key]["counts"] += intersect
            bed_dict[bed_key]["window_counts"] += 1
        else:
            dict_in={"length":length,"counts":intersect,"window_counts":1.0}
            bed_dict[bed_key] = dict_in
      
    return bed_dict


def convert_dictToBedIntersected(target_dict, length=0, counts=1, window_counts=2):


    for target_keys, target_values in target_dict.items():
        
        temp = []
        info = target_keys.split("_cut_")

        temp.append(info[1])
        temp.append("0")
        temp.append(str(list(target_values.values())[length]))
        temp.append(info[0])
        temp.append(str(info[2]))
        temp.append(info[3])
        temp.append(str(list(target_values.values())[counts]/list(target_values.values())[window_counts]))

        #print(temp)
        out.write("\t".join(temp) + "\n")
    
    out.close()
    return


step1 = bedIntersected_dictionary()
convert_dictToBedIntersected(step1)

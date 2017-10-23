1

import re
import math
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Arguments for function..')

parser.add_argument('-i', dest='infile', type=str)
parser.add_argument('-c', dest='chains', type=str, nargs='+')
parser.add_argument('-o', dest='outfile', type=str)
parser.add_argument('-r', dest='reference', type=str)
parser.add_argument('-co', dest='cutoff', type=float)
args = parser.parse_args()

#def remove_chains(target_file, target_chains, output_file):
#    with open(target_file, 'r') as ifile:
#        for line in ifile:
#            els = line.split()
#            for i in target_chains:
#                if els[0] == 'ATOM' or els[0] == 'HETATM':
#                    if re.match(i, els[4]) or re.match(i + '/d', els[4]):
#                        with open(output_file, 'a') as ofile:
#                            ofile.write(line + "\n")

def water_compare(target_file, reference_file, cutoff):

    with open(target_file, 'r') as ofile:
        conv_list = []
        ref_list = []
        dist_list = []
        for line in ofile:
            if re.match('(.*)HOH(.*)', line):
                line = line.split()
                X1 = float(line[6])
                Y1 = float(line[7])
                Z1 = float(line[8])
                Conserved = int(line[5])

                with open(reference_file, 'r') as rfile:
                    for line in rfile:
                        if re.match('(.*)HOH(.*)', line):
                            line = line.split()
                            X2 = float(line[6])
                            Y2 = float(line[7])
                            Z2 = float(line[8])
                            pot_ref = int(line[5])

                            distance =float(math.sqrt(((X1-X2) ** 2) + ((Y1-Y2) ** 2) + ((Z1-Z2) ** 2)))
                            

                            if distance <= cutoff:
                                conv_list.append(Conserved)
                                ref_list.append(pot_ref)
                                dist_list.append(distance)
                                
    return conv_list, ref_list, dist_list 

#TODO: Find a pure python way to do this
def remove_duplicates(clist, dlist, rlist):
    data = (zip(clist, rlist, dlist))
    df = pd.DataFrame(data,columns=["Conserved ID", "Reference ID", "Distance"])
    df.sort_values("Distance", inplace=True, ascending=False)
    df.drop_duplicates(subset=["Conserved ID"], keep="last", inplace=True)
    return df["Conserved ID"].tolist(), df["Reference ID"].tolist(), df["Distance"].tolist()

def generate_ref_pdb(reference_file, pred_list, output_file, input_file):
    ofile = open(output_file, 'w')
    ofile.seek(0)
    with open(input_file) as ifile:
        for line in ifile:
            if re.match('(.*)HOH(.*)', line):
                for i in pred_list:
                    if int(line.split()[5]) == int(i):
                        ofile.write(line + '\n')
            else:
                ofile.write(line + '\n')
                
    ofile.close()
    
clist, rlist, dlist = water_compare(args.infile, args.reference, args.cutoff)
dd_clist, dd_rlist, dd_dlist = remove_duplicates(clist, dlist, rlist)
generate_ref_pdb(args.reference, dd_rlist, args.outfile, args.infile)

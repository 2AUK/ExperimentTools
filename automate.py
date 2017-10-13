"""
Automation of Predicted vs Experimental Comparison
|Requires a specific layout of the filesystem before using
|otherwise it won't work. Main reason for this is because it's easier to work
|with this setup
TODO: Generalise this. Don't need a specific filesystem layout, just the pdb
files and their paths
"""

import os
from math import sqrt
import re
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter

DEST = "/Users/AbdullahAhmad/Desktop/Aspartic_Proteases_Automation"
FOLDERS = os.walk(DEST).next()[1]
OUTPUT = "/Users/AbdullahAhmad/Desktop/Aspartic_Proteases_Automation_Output"
DIST_CUTOFF = 2
DENS_CUTOFF = 5
PROTEIN_COUNT = len(FOLDERS)
REFERENCE_STRUCT = '4CMS'

def plot_func(data_list, ref):
    """
    THIS IS SO MESSY
    """
    combined = []
    reference = []
    for main, comp, percentage in zip(*[iter(data_list)]*3):
        if comp == "Combined":
            combined.append(percentage)
        else:
            reference.append(percentage)
    ind = np.arange(PROTEIN_COUNT)
    width = 0.25
    fig, ax = plt.subplots(figsize=(20,10))
    rects1 = ax.bar(ind, combined, width, color='r')
    rects2 = ax.bar(ind+width, reference, width, color='y')
    ax.set_ylabel('Percentage Correct with ' + str(DENS_CUTOFF) + ' Cut-off')
    ax.set_title(ref)
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(FOLDERS)
    ax.legend((rects1[0], rects2[0]), ('Combined', 'Reference: ' + ref))
    plt.show()
    fig.savefig("Score" + ref + ".png")
    
        
def compare(pids, tdir, odir, cutoff, ref):
    """
    Main Comparison function
    STARTING REWRITE MY FRIEND
    """
    tot_list = []
    for pid in pids:
        current_ref = (tdir + "/" + pid + "/" + pid + "_CW.pdb")
        if not os.path.exists(odir + "/" + pid):
            os.mkdir(odir + "/" + pid)
        for auto in pids:
            conserved = []
            predicted = []
            density = []
            distance = []
            generated_file = (tdir + "/" + auto + "/" + "6-Comparison/O.pdb")
            with open(current_ref, "r") as exp_file:
                for line in exp_file:
                    if re.match('(.*)HOH(.*)', line):
                        els = line.split()
                        x_1 = float(els[6])
                        y_1 = float(els[7])
                        z_1 = float(els[8])
                        expid = int(els[5])
                        with open(generated_file, 'r') as pred_file:
                            for line in pred_file:
                                if re.match('(.*)HETATM(.*)', line):
                                    subels = line.split()
                                    x_2 = float(subels[5])
                                    y_2 = float(subels[6])
                                    z_2 = float(subels[7])
                                    predid = int(subels[1])
                                    dens = float(subels[8])
                                    dist = float(sqrt((float(x_1 - x_2) ** 2) + \
                                    (float(y_1 - y_2) ** 2) + \
                                    (float(z_1 - z_2) ** 2)))
                                    if dist <= cutoff:
                                        conserved.append(expid)
                                        predicted.append(predid)
                                        density.append(dens)
                                        distance.append(dist)
                                        elist, plist, denlist, dlist = remove_duplicates(conserved, predicted, density,\
                                                distance, odir, pid, auto)
            a = (float(len(plist)) / float(max(plist))) * 100
            if a != 0 and ("Combined" in pred_file.name or ref in pred_file.name): 
                tot_list.append(pid), tot_list.append(auto), tot_list.append(a)
    return tot_list

def init_data(input_path):
    """
    Initialise a dictionary of protein IDs and their directories for reading
    |Arguments
    |input_path - a directory containing the relevant proteins. Currently going 
    |            to use the numbered-folder setup.
    """
    #Scan folders in the input directory
    protein_ids = os.walk(input_path).next()[1]
    for prot_id in protein_ids: pass
        

def vector_distance_total(input_list):
    """
    Calculating distance, using a generator
    |Arguments
    |input_list - an list of floats where every 6 values correspond to 2 vectors
    """
    for x_1, y_1, z_1, x_2, y_2, z_2 in zip(*[iter(input_list)]*6):
        x_m_x, y_m_y, z_m_z = x_1 - x_2, y_1 - y_2, z_1 - z_2
        distance = sqrt(x_m_x ** 2 + y_m_y ** 2 + z_m_z ** 2) 
        yield distance
        
def sortlist(input_list, col, direc):
    """
    Function for sorting rows in a table
    |Arguments
    |input_list - list of lists, make sure the elements of the sub lists are numbers
    |             not strings
    |col - which column to sort in the table
    |direc - ascending, descending
    """
    return sorted(input_list, key=lambda row: (row[col]), reverse=direc)

        
def duplicate_filter(input_list):
    """
    Pure python method to remove duplicates from a table of data.
    |Arguments
    |input_list - a list of lists, with each list being a row in a table
    """
    sortedlist = sortlist(input_list, 3, False)
    outputlist = []
    seen = set()
    for row in sortedlist:
        if tuple(row[0:3]) in seen:
            continue
        seen.add(tuple(row[0:3]))
        outputlist.append(row)
    for i in sortlist(outputlist, 1, False):
        yield i
        
def process_data(proteins, output):
    """
    Read relevant data from each protein and process. This is the main function.
    |Arguments
    |proteins - a dictionary of the PDB ID and filesystem location
    |output - some path to output all the data to, either as .txt or maybe .csv
    """
    pass
                            





f = open("/Users/AbdullahAhmad/Desktop/Aspartic_Proteases_Automation_Output/1AM5/Ref_1AM5_CF_Combined.txt")
f.readline()
inputs = []
for line in f:
    line = map(float, line.split())
    line[0] = int(line[0])
    line[1] = int(line[1])
    inputs.append(line)
                 
new = duplicate_filter(inputs)

for i in new:
    print "\t".join(map(str, i))
    